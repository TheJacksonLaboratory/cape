#' Runs marker regression on each individual genetic marker
#' 
#' This function performs marker regression to associate
#' individual markers with traits (or eigentraits).
#' If n_perm is greater than 0, permutations are run to 
#' determine effect size thresholds for the alpha values
#' provided. The default alpha values are 0.05 and 0.01.
#' Covariates are specified in the cape parameter file.
#' 
#' model_family indicates the model family of the phenotypes
#' This can be either "gaussian" or "binomial". If this argument
#' is length 1, all phenotypes will be assigned to the same
#' family. Phenotypes can be assigned different model families by
#' providing a vector of the same length as the number of phenotypes,
#' indicating how each phenotype should be modeled.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object.
#' @param kin_obj a kinship object. If NULL, the kinship correction is not performed.
#' @param n_perm integer number of permutations. Permutation results are only
#' used in \code{\link{plot_singlescan}}. They are not used for any other piece
#' of the Cape analysis and may be safely omitted. The default number of permutations
#' is 0.
#' @param alpha significance level if permutations are being run. If permutations are
#' run effect size thresholds for each alpha level are cacluated using the extreme value
#' distribution.
#' @param model_family A vector indicating the model family of the phenotypes. This can 
#' be either "gaussian" or "binomial." If length 1, all phenotypes will be assigned to the 
#' same family. Phenotypes can be assigned different model families by
#' providing a vector of the same length as the number of phenotypes,
#' indicating how each phenotype should be modeled.
#' @param run_parallel Whether to run on parallel CPUs
#' @param n_cores The number of CPUs to use if run_parallel is TRUE
#' @param verbose Whether to print progress to the screen
#' @param overwrite_alert Used 
#'
#' @return Returns a list of the singlescan results. The list is
#' of length seven, and has the following elements: 
#'    alpha: The alpha values set in the argument alpha
#'    alpha_thresh: The calculated effect size thresholds at each alpha if permutations are run.
#'    ref_allele: The allele used as the reference allele
#'    singlescan_effects: The effect sizes (beta coefficients) from the single-locus linear models
#'    singlescan_t_stats: The t statistics from the single-locus linear models
#'    locus.p_vals: Marker-level p values
#'    locus_score_scores: Marker-level test statistics.
#'
#' @seealso \code{\link{plot_singlescan}}
#' 
#' @import parallel
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom Matrix rankMatrix
#' @importFrom stats var
#' @importFrom stringr str_replace
#' 
#' @examples 
#' \dontrun{
#' singlescan_obj <- singlescan(data_obj, geno_obj, kin_obj)
#' }
#' 
#' @export
#' 
singlescan <- function(data_obj, geno_obj, kin_obj = NULL, n_perm = 0, 
  alpha = c(0.01, 0.05), model_family = "gaussian", run_parallel = FALSE, 
  n_cores = 4, verbose = FALSE, overwrite_alert = TRUE) {
  
  ref_allele <- data_obj$ref_allele
  scan_what <- data_obj$scan_what
  #use_kinship <- data_obj$use_kinship
  use_kinship <- as.logical(length(kin_obj))

  if(!run_parallel){n_cores = 1}
  

  if(overwrite_alert){
    choice <- readline(prompt = "Please make sure you assign the output of this function to a singlescan_obj, and NOT the data_obj. It will overwrite the data_obj.\nDo you want to continue (y/n) ")
    if(choice == "n"){stop()}
  }
  
  if(length(model_family) == 1){
    model_family <- rep(model_family, ncol (data_obj$pheno))
  }
  
  data_obj$model_family <- model_family
  
  #===============================================================
  gene <- get_geno(data_obj, geno_obj)
  pheno <- get_pheno(data_obj, scan_what)	
  n_phe = dim(pheno)[2]
  chr_which <- unique(data_obj$chromosome)
  
  #get the covariates and assign the variables
  #to the local environment
  covar_info <- get_covar(data_obj)
  covar_names <- covar_info$covar_names
  covar_table <- covar_info$covar_table
  for(i in 1:length(covar_info)){
    assign(names(covar_info)[i], covar_info[[i]])
  }
  
  n_covar <- length(covar_names)
  
  #perform a check of covariates
  if(!is.null(covar_table)){
    cov_var <- apply(covar_table, 2, var)
    zero_locale <- which(cov_var == 0)
    if(length(zero_locale) > 0){
      stop(paste(covar_names[zero_locale], "has zero variance."))
    }		
    
    # TODO this code is duplicated in get_eigentraits
    
    #also remove the NAs and check the matrix for rank
    not_na_locale <- which(!is.na(rowSums(covar_table)))
    no_na_cov <- as.array(covar_table[not_na_locale,,drop=FALSE])
    design_cov <- cbind(rep(1, dim(no_na_cov)[1]), no_na_cov)
    rank_cov <- rankMatrix(design_cov)
    if(rank_cov[[1]] < dim(design_cov)[2]){
      stop("The covariate matrix does not appear to be linearly independent.\nIf you are using dummy variables for groups, leave one of the groups out.")
    }
  }
  
  if(is.null(n_perm) || n_perm < 2){alpha = "none"}
  #===============================================================
  
  singlescan_obj <- vector(mode = "list", length = 7)
  names(singlescan_obj) <- c("alpha", "alpha_thresh", "ref_allele", "singlescan_effects", "singlescan_t_stats", "locus.p_vals", "locus_score_scores")
  
  #===============================================================
  singlescan_obj$covar <- covar_names
  singlescan_obj$alpha <- alpha
  #===============================================================
  
  
  #==================================================================
  #if we are using a kinship correction, make sure the phenotypes
  #are mean-centered and there are no missing values in the genotype
  #matrix.
  #==================================================================
  if(!is.null(kin_obj)){
    pheno_means <- apply(pheno, 2, mean)
    tol = 0.01
    non_zero <- intersect(which(pheno_means > 0+tol), which(pheno_means < 0-tol))
    if(length(non_zero) > 0){
      warning("Phenotypes must be mean-centered before performing kinship corrections.")
      message("Mean-centering phenotypes using norm_pheno()")
      data_obj <- norm_pheno(data_obj)
    }
    
   # missing_vals <- which(is.na(gene))
   # if(length(missing_vals) > 0){
   #   warning("There are missing values in the genotype matrix. Please use impute_missing_geno().")
   #   data_obj <- impute_missing_geno(data_obj, geno_obj = geno_obj, run_parallel = run_parallel, n_cores = n_cores)["data_obj"]
   # }
  }
  #==================================================================
  
  #Get the dimension names to minimize confusion	
  geno_dims <- get_geno_dim()
  mouse_dim <- geno_dims[which(names(geno_dims) == "mouse")]
  allele_dim <- geno_dims[which(names(geno_dims) == "allele")]
  locus_dim <- geno_dims[which(names(geno_dims) == "locus")]
  
  n_phe <- dim(pheno)[2]
  
  #first do the permutations to get the significance threshold
  #results will be compared to the significance threshold to 
  #determine which markers to use as covariates
  
  if(n_perm > 0){
    if(verbose){cat("\n\nPerforming permutations to calculate significance threshold...\n")}			
    singlescan_obj$alpha_thresh <- genome_wide_threshold_1D(data_obj, geno_obj, 
    n_perm = n_perm, scan_what = scan_what, ref_allele = ref_allele, alpha = alpha, 
    model_family = model_family, n_cores = n_cores, run_parallel = run_parallel,
    verbose = verbose)
  }else{
    if(verbose){cat("Not performing permutations in singlescan.\n")}
  }
  
  
  #check for a reference allele, pull it out of the 
  #allele names here and add it to the data object
  if(length(ref_allele) != 1){ #add a check for the reference allele
    stop("You must specify one reference allele")
  }
  ref_col <- which(dimnames(gene)[[allele_dim]] == ref_allele)
  if(length(ref_col) == 0){
    stop("I can't find reference allele: ", ref_allele)
  }
  new_allele_names <- dimnames(gene)[[allele_dim]][-ref_col]
  singlescan_obj$ref_allele <- ref_allele	
  
  #=====================================================================
  #internal functions
  #=====================================================================
  
  #This function takes the results from get_stats_multiallele
  #and parses them into the final arrays
  add_results_to_array <- function(result_array, results_list, stat_name, is_covar = FALSE){
    row_num <- which(rownames(results_list[[1]][[1]]) == stat_name)
    #find the next spot to put data
    #Each successive phenotype is stored
    #in the 2nd dimension of the array
    next_spot_locale <- min(which(is.na(result_array[nrow(result_array),,1])))
    result_mat <- t(sapply(results_list, function(x) as.vector(x[[1]][row_num,])))
    
    if(is_covar){ #if we are looking at a covariate, we need to expand it
      result_mat <- matrix(result_mat, nrow = ncol(result_mat), ncol = dim(result_array)[3])
    }
    
    placement_locale <- match(names(results_list), rownames(result_array))
    # na_locale <- which(is.na(result_array[,next_spot_locale,1]))
    result_array[placement_locale,next_spot_locale,] <- result_mat
    return(result_array)	
  }
  
  add_flat_results_to_array <- function(result_array, model, pheno_num){
    betas <- as.vector(model[[2]]$beta)
    num_markers <- dim(result_array)[1]
    num_alleles <- dim(result_array)[3]
    start_pos <- 1
    for(i in 1:num_markers){
      locus_results <- betas[start_pos:(start_pos+num_alleles-1)]
      result_array[i,pheno_num,] <- locus_results
      start_pos = start_pos + num_alleles
    }
    return(result_array)
  }
  
  
  #=====================================================================
  #end of internal functions
  #=====================================================================
  
  
  #=====================================================================
  #begin code for multi-allelic cross
  #=====================================================================
  #In the multi-allele case, we want to collect
  #three 3D arrays each of num_marker by num_pheno by num_allele:
  #array of t statistics (for plotting p vals of regressions)
  #array of effects (betas) (for effect plots)
  #array of covar flags (for use in pair.scan)
  
  t_stat_array <- array(dim = c(dim(gene)[[locus_dim]]+n_covar, dim(pheno)[2], (dim(gene)[[allele_dim]]-1)))
  effect_array <- array(dim = c(dim(gene)[[locus_dim]]+n_covar, dim(pheno)[2], (dim(gene)[[allele_dim]]-1)))
  dimnames(t_stat_array) <- dimnames(effect_array) <- list(c(dimnames(gene)[[locus_dim]], covar_names), dimnames(pheno)[[2]], new_allele_names)
  
  #make arrays to hold locus-by-locus stats
  locus_score_scores <- matrix(NA, ncol = n_phe, nrow = dim(gene)[[locus_dim]]+n_covar)
  locus_p_values <- matrix(NA, ncol = n_phe, nrow = dim(gene)[[locus_dim]]+n_covar)
  rownames(locus_score_scores) <- rownames(locus_p_values) <- c(dimnames(gene)[[locus_dim]], covar_names)
  colnames(locus_score_scores) <- colnames(locus_p_values) <- colnames(pheno)
  for(i in 1:n_phe){
    if(verbose){cat("\nScanning trait:", colnames(pheno)[i], "\n")}
    #take out the response variable
    phenotype <- pheno[,i,drop=FALSE]
    ph_family = model_family[i]
    #get corrected genotype and phenotype values for each phenotype-chromosome pair
    if(use_kinship){
      sink(file.path(data_obj$results_path,"regress.warnings.singlescan")) #create a temporary output file for the regress warnings
      # TODO check if dim(kin.obj)[1] == length(phenoV) == length(covarV) when using covariates
      if(data_obj$kinship_type == "ltco"){
        cor_data <- lapply(chr_which, function(x) kin_adjust(kin_obj, gene, 
        chr1 = x, chr2 = x, phenoV = phenotype, covarV = covar_table))
        names(cor_data) <- chr_which
      }else{
        cor_data <- vector(mode = "list", length = 1)
        cor_data[[1]] <- kin_adjust(kin_obj, gene, phenoV = phenotype, 
        covarV = covar_table)
      }
      sink(NULL)
    }else{
      cor_data <- vector(mode = "list", length = 1)
      cor_data[[1]] <- list("corrected_pheno" = phenotype, "corrected_geno" = gene, 
      "corrected_covar" = covar_table)
    }
    
    results_by_chr <- vector(mode = "list", length = length(cor_data))							
    
    for(ch in 1:length(cor_data)){
      if(use_kinship && verbose){cat(" Chr", ch, "... ", sep = "")}
      if(length(cor_data) == 1){chr_locale <- 1:dim(gene)[3]}
      if(length(cor_data) > 1){chr_locale <- which(data_obj$chromosome == names(cor_data)[ch])}
      c_geno <- cor_data[[ch]]$corrected_geno[,,chr_locale,drop=FALSE]
      c_pheno <- cor_data[[ch]]$corrected_pheno
      c_covar <- cor_data[[ch]]$corrected_covar
      
      if (run_parallel) {
        cl <- makeCluster(n_cores)
        registerDoParallel(cl)
        
        # the following line adds package variables to the parallel worker environments
        # copy functions in the package to the workers
        # TODO remove this hardcoded line, supply a variable to the Cape.obj containing the full path
        #cape_dir <- "/Users/ramamg/Desktop/JAX/Projects/CAPE/cape/cape_pkg"
        cape_dir_full <- find.package("cape")
        cape_dir <- str_replace(cape_dir_full,"cape_pkg/cape","cape_pkg")
        clusterExport(cl, "cape_dir", envir=environment())
        clusterEvalQ(cl, .libPaths(cape_dir))
        results_by_chr <- foreach(x = 1:dim(c_geno)[locus_dim], .packages = 'cape') %dopar% {
          # Note that "Show Diagnostics" in RStudio will throw a warning that the `x` variable below is undefined
          # but it actually is defined in the foreach line above. You can safely ignore the warning.
          cape::get_stats_multiallele(phenotype = c_pheno, genotype = c_geno[,,x], covar_table = c_covar, ph_family, ref_col)
        }
        stopCluster(cl)
        
      } else {
        
        results_by_chr <- c()
        index <- 1:dim(c_geno)[locus_dim]
        for (x in index) {
          results_by_chr[[x]] <- get_stats_multiallele(phenotype = c_pheno, 
          genotype = c_geno[,,x], covar_table = c_covar, ph_family, ref_col)
        }
        
      }
      names(results_by_chr) <- dimnames(c_geno)[[locus_dim]]	
      
      t_stat_array <- add_results_to_array(result_array = t_stat_array, results_list = results_by_chr, stat_name = "t_stat")
      effect_array <- add_results_to_array(result_array = effect_array, results_list = results_by_chr, stat_name = "slope")
      locus_score_scores[chr_locale,i] <- unlist(lapply(results_by_chr, function(x) x$score))
      
    } #end looping through data corrected by chromosome (loco)
    
    #if there are covariates, run them through too
    if(!is.null(covar_table)){
      if(verbose){cat("\nTesting covariates \n")}
      #get corrected genotype and phenotype values for the overall kinship matrix
      if(use_kinship){
        sink(file.path(data_obj$results_path,"regress.warnings.covar"))
        # TODO check if dim(kin.obj)[1] == length(phenoV) == length(covarV) when using covariates
        cor_data <- kin_adjust(kin_obj, gene, phenoV = phenotype, covarV = covar_table, 
        verbose = verbose)
        sink(NULL) #stop sinking to the file
      }else{
        cor_data <- list("corrected_pheno" = phenotype, "corrected_geno" = gene, "corrected_covar" = covar_table)
      }
      c_geno <- cor_data$corrected_geno
      c_pheno <- cor_data$corrected_pheno
      c_covar <- cor_data$corrected_covar 
      covar_results <- apply(c_covar, 2, function(x) get_stats_multiallele(c_pheno, x, c_covar, ph_family, ref_col))
      names(covar_results) <- data_obj$p_covar
      t_stat_array <- add_results_to_array(result_array = t_stat_array, results_list = covar_results, stat_name = "t_stat", is_covar = TRUE)
      effect_array <- add_results_to_array(effect_array, covar_results, "slope", is_covar = TRUE)
      first_na <- min(which(is.na(locus_score_scores[,i])))
      locus_score_scores[(first_na):(first_na+n_covar-1),i] <- unlist(lapply(covar_results, function(x) x$score))
    } #end case for an existing covariate table
    
  }      #end looping through phenotypes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  if(verbose){cat("\n")}
  
  singlescan_obj$singlescan_effects <- effect_array
  singlescan_obj$singlescan_t_stats <- t_stat_array
  # singlescan_obj$locus.p_vals <- locus_p_values
  singlescan_obj$locus_score_scores <- locus_score_scores
  singlescan_obj$n_perm <- n_perm
  
  #unlink(file.path(data_obj$results_path,"regress.warnings")) #remove the temporary file
  return(singlescan_obj)
  
}

