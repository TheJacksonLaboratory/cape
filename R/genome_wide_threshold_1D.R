#' Calculate a genome-wide significance threshold for the single-variant scan
#' 
#' This function uses permutation testing to calculate a genome-wide significance
#' threshold for the single-variant scan. Two user-defined thresholds are 
#' calculated: the lower threshold (higher alpha) can be used to determine which 
#' variants will be used in the pairwise scan. The higher threshold (lower alpha) 
#' is used to determine which variants are used as covariates in the pairwise scan. 
#' In each permutation, the phenotype or eigentrait is shuffled and all markers 
#' are retested with the permuted phenotype. The regression coefficients are 
#' collected from each permutation and the extreme value distribution is used to 
#' determine thresholds for the user-defined alpha values.
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object.
#' @param n_perm The number of permutations to perform. The default is 100.
#' @param scan_what A character value that uniquely specifies whether the 
#' eigentraits or phenotypes should be scanned. Options are "eigentraits"
#' or "raw_traits".
#' @param ref_allele \code{\link{Cape}} requires that one of the alleles 
#' in the population be selected as a reference allele. The effects of all 
#' alleles are then reported as the effects relative to the reference allele. 
#' In the DO population B represents the C57BL/6J (B6) mouse. Because this strain is 
#' considered the standard strain, the default reference allele is B. 
#' If another allele has a positive effect, it means that it increases the given 
#' phenotype relative to that in the B6 mouse. Any other allele can be selected 
#' as the reference allele simply by specifying it with this argument. 
#' @param alpha The alpha value(s) used to calculate significance levels. This should
#' be a vector of any length of numerical values between 0 and 1. The default is
#' a vector of length two: c(0.01, 0.05).
#' @param model_family Indicates the model family of the phenotypes. This can be 
#' either "gaussian" or "binomial".
#' @param run_parallel A logical value indicating whether the process should
#' be run in parallel. Defaults to FALSE.
#' @param n_cores integer number of cores to use if running in parallel
#' @param verbose A logical value indicating whether to print progress to the screen.
#' Defaults to FALSE.
#' 
#' @return Returns a vector the same length as alpha indicating the
#' thresholds for each value of alpha.
#' 
#' @import parallel
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom stats coef glm
#' @importFrom evd qgev fgev
#' @keywords internal
#' 
genome_wide_threshold_1D <- function(data_obj, geno_obj, n_perm = 100, 
                                     scan_what = c("eigentraits", "raw_traits"), 
                                     ref_allele = NULL, alpha = c(0.01, 0.05), 
                                     model_family, run_parallel = FALSE, n_cores = 4,
                                     verbose = FALSE){
  
  if(!run_parallel){n_cores = 1}
    
  if(n_perm < 2){
    stop("You must do more than one permutation.")
  }
  
  gene <- get_geno(data_obj, geno_obj)
  
  #Get the dimension names to minimize confusion
  geno_dims <- get_geno_dim()
  mouse_dim <- geno_dims[which(names(geno_dims) == "mouse")]
  allele_dim <- geno_dims[which(names(geno_dims) == "allele")]
  locus_dim <- geno_dims[which(names(geno_dims) == "locus")]
   
  #check to make sure a reference allele has been chosen
  #and that we can find it in the array
  if(length(ref_allele) != 1){
    stop("You must pick one reference allele")
  }else{
    ref_col <- which(dimnames(gene)[[allele_dim]] == ref_allele)
    if(length(ref_col) < 1){
      stop("I couldn't find reference allele: ", ref_allele)
    }
  }
  
  
  #calculate the numbers of markers, phenotypes and samples
  # TODO remove the lines below if not needed
  # n_gene <- dim(gene)[locus_dim]
  # n_allele <- dim(gene)[allele_dim]-1 #find the number of alleles we will be regressing on
  
  #pull out genotype and phenotype data for
  #clarity of code later.
  #If the user has not specified a scan_what,
  #from the outer function (singlescan.R),
  #default to eigentraits, basically, if eigen,
  #et, or ET are anywhere in the string, use the
  #eigentrais, otherwise, use raw phenotypes
  use_eigentraits <- length(c(grep("eigen", scan_what), grep("ET", scan_what), grep("et", scan_what)))
  if(use_eigentraits){
    pheno <- data_obj$ET
  }else{
    pheno <- data_obj$pheno
  }
  
  
  num_samples <- dim(pheno)[1]
  
  if(length(model_family) == 1){
    model_family <- rep(model_family, ncol(pheno))
  }
  
  
  #============================================================	
  #functions
  #============================================================	
  #do the regression at a single locus and return the model
  multi_allele_regress <- function(locus_mat, phenotype, model_family){
    #remove the reference allele from the matrix
    regress_mat <- locus_mat[,-ref_col]
    model <- glm(phenotype~regress_mat, family = model_family)
    model_coef <- coef(summary(model))
    if(nrow(model_coef) == 2){
	    #gather the t statistics for all alleles
    		t_stats <- as.vector(model_coef[2:dim(model_coef)[1], 3])
    		}else{
    		t_stats <- NA
    		}
    return(t_stats)
  }
  
  get_s <- function(evd_result, alpha){
    s <- qgev(1-alpha,loc=evd_result$estimate[1], scale=evd_result$estimate[2], shape=evd_result$estimate[3], lower.tail = TRUE)
    return(s)
  }
  
  
  #This funtion runs a single permutation
  #it can be parallelized to speed things up
  one_perm <- function(){
    #shuffle the vector of individuals
    #This keeps all correlation structure
    #between phenotypes and genotypes the
    #same, but reassigns which phenotypes
    #are associated with which genotypes
    sampled_vector <- sample(1:num_samples)
    
    #permute the individuals
    gene_perm <- gene[sampled_vector,,]
    
    #for each phenotype, do the regression on all loci
    #at once and pull out the t statistic
    
    stat_mat <- array(NA, dim = c(dim(gene_perm)[locus_dim], ncol(gene_perm)-1, ncol(pheno)))
    
    # ptm <- proc.time()
    for(et in 1:ncol(pheno)){	 	  			
      stat_mat[,,et] <- Reduce("rbind", lapply(1:dim(gene_perm)[locus_dim], function(x) multi_allele_regress(locus_mat = gene_perm[,,x], phenotype = pheno[,et], model_family = model_family[et])))
    }
    
    #find the maximum t statistic for the permutation
    #and record it.
    max_stat <- apply(stat_mat, 2, function(x) max(x, na.rm = TRUE))
    return(max_stat)        				
  }
  
  #====================================================
  
  # TODO remove the two lines below if not needed
  #create a matrix to hold permutation results
  #perm_max <- array(NA, dim = c(n_perm, ncol(gene)-1, ncol(pheno)))
  
  
  if (run_parallel) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    cape_dir_full <- find.package("cape")
    cape_dir <- str_replace(cape_dir_full,"cape_pkg/cape","cape_pkg")
    clusterExport(cl, "cape_dir", envir=environment())
    clusterEvalQ(cl, .libPaths(cape_dir))
    max_stat <- foreach(p = 1:n_perm, .combine = "rbind")  %dopar% {
      one_perm()
    }
    stopCluster(cl)
    
  } else {
    
    max_stat <- c()
    index <- 1:n_perm
    for (p in index) {
      if(verbose){report_progress(p, n_perm)}
      max_stat <- rbind(max_stat, one_perm())
    }
    
  }
  # proc.time() - ptm
  
  #====================================================
  # calculate the pairscan and covar thresholds and 
  # add them to the data object.
  #====================================================
  
  #apply the extreme value distribution to the results
  evd <- apply(max_stat, 2, function(x) fgev(x, std.err = FALSE))
  data_obj$save_rds(max_stat, "singlescan_permutations.RDS")
  
  s <- vector(mode = "list", length = length(alpha))
  for(a in 1:length(alpha)){
    s[[a]] <- as.vector(sapply(evd, function(x) get_s(x, alpha[a])))
  }
  
  #calculate one threshold over all phenotypes
  thresholds <- lapply(s, mean)
  names(thresholds) <- alpha
  
  if(verbose){cat("\n")}

  return(thresholds)
  
  
}
