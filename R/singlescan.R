#' Runs a 1D scan of the phenotype data to assess which markers to use as covariates for each phenotype
#' 
#' This script performs marker regression to associate
#' individual markers with either traits or eigentraits.
#' If n.perm is greater than 0, permutations are run to 
#' determine effect size thresholds for the alpha values
#' provided. The default alpha values are 0.05 and 0.01.
#' Covariates are specified in the cape parameter file.
#' 
#' model.family indicates the model family of the phenotypes
#' This can be either "gaussian" or "binomial"
#' if length 1, all phenotypes will be assigned to the same
#' family. If the phenotypes are a mixture, model.family can
#' be a vector of length Np, where Np is the number of phenotypes
#' indicating which phenotype belongs in which family of distributions.
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object. If this is not supplied then it is generated here.
#' @param kin.obj a kinship object
#' @param n.perm integer number of permutations (default=100)
#' @param alpha real significance level (either 0.01 or 0.05)
#' @param model.family array indicating the model family of the phenotypes. This can be either "gaussian" or "binomial" if length 1, all phenotypes will be assigned to the same family. If the phenotypes are a mixture, model.family can be a vector of length Np, where Np is the number of phenotypes indicating which phenotype belongs in which family.
#' @param run.parallel boolean, default = TRUE
#' @param n.cores integer number of cores to use if running in parallel
#' @param verbose boolean
#' @param overwrite.alert boolean, for use when the script is run from the command line
#'
#' @return \code{list("data.obj" = data.obj, "geno.obj" = geno.obj)}
#'
#' @export
singlescan <- function(data.obj, geno.obj, kin.obj = NULL, n.perm = 100, alpha = c(0.01, 0.05), 
                       model.family = "gaussian", run.parallel = TRUE, n.cores = 4, verbose = FALSE, 
                       overwrite.alert = TRUE) {
  
  ref.allele <- data.obj$ref_allele
  scan.what <- data.obj$scan_what
  
  if(!run.parallel){n.cores = 1}
  
  if(!is.null(kin.obj)){use.kinship = TRUE}
  if(is.null(kin.obj)){use.kinship = FALSE}
  
  if(overwrite.alert){
    choice <- readline(prompt = "Please make sure you assign the output of this function to a singlescan.obj, and NOT the data.obj. It will overwrite the data.obj.\nDo you want to continue (y/n) ")
    if(choice == "n"){stop()}
  }
  
  if(length(model.family) == 1){
    model.family <- rep(model.family, ncol (data.obj$pheno))
  }
  
  data.obj$model_family <- model.family
  
  #===============================================================
  gene <- get.geno(data.obj, geno.obj)
  pheno <- get.pheno(data.obj)	
  n.phe = dim(pheno)[2]
  chr.which <- unique(data.obj$chromosome)
  
  #get the covariates and assign the variables
  #to the local environment
  covar.info <- get.covar(data.obj)
  covar.names <- covar.info$covar.names
  covar.table <- covar.info$covar.table
  for(i in 1:length(covar.info)){
    assign(names(covar.info)[i], covar.info[[i]])
  }
  
  n.covar <- length(covar.names)
  
  #perform a check of covariates
  if(!is.null(covar.table)){
    cov.var <- apply(covar.table, 2, var)
    zero.locale <- which(cov.var == 0)
    if(length(zero.locale) > 0){
      stop(paste(covar.names[zero.locale], "has zero variance."))
    }		
    
    # TODO this code is duplicated in get.eigentraits
    
    #also remove the NAs and check the matrix for rank
    not.na.locale <- which(!is.na(rowSums(covar.table)))
    no.na.cov <- as.array(covar.table[not.na.locale,,drop=FALSE])
    design.cov <- cbind(rep(1, dim(no.na.cov)[1]), no.na.cov)
    rank.cov <- Matrix::rankMatrix(design.cov)
    if(rank.cov[[1]] < dim(design.cov)[2]){
      stop("The covariate matrix does not appear to be linearly independent.\nIf you are using dummy variables for groups, leave one of the groups out.")
    }
  }
  
  if(is.null(n.perm) || n.perm < 2){alpha = "none"}
  #===============================================================
  
  singlescan.obj <- vector(mode = "list", length = 7)
  names(singlescan.obj) <- c("alpha", "alpha.thresh", "ref.allele", "singlescan.effects", "singlescan.t.stats", "locus.p.vals", "locus.score.scores")
  
  #===============================================================
  singlescan.obj$covar <- covar.names
  singlescan.obj$alpha <- alpha
  #===============================================================
  
  
  #==================================================================
  #if we are using a kinship correction, make sure the phenotypes
  #are mean-centered and there are no missing values in the genotype
  #matrix.
  #==================================================================
  if(!is.null(kin.obj)){
    pheno.means <- apply(pheno, 2, mean)
    tol = 0.01
    non.zero <- intersect(which(pheno.means > 0+tol), which(pheno.means < 0-tol))
    if(length(non.zero) > 0){
      warning("Phenotypes must be mean-centered before performing kinship corrections.")
      cat("Mean-centering phenotypes using norm.pheno()")
      data.obj <- norm.pheno(data.obj)
    }
    
    missing.vals <- which(is.na(gene))
    if(length(missing.vals) > 0){
      warning("There are missing values in the genotype matrix. Please use impute.missing.geno().")
      data.obj <- impute.missing.geno(data.obj, geno.obj)["data.obj"]
    }
  }
  #==================================================================
  
  #Get the dimension names to minimize confusion	
  mouse.dim <- which(names(dimnames(gene)) == "mouse")
  locus.dim <- which(names(dimnames(gene)) == "locus")
  allele.dim <- which(names(dimnames(gene)) == "allele")
  
  n.phe <- dim(pheno)[2]
  
  #first do the permutations to get the significance threshold
  #results will be compared to the significance threshold to 
  #determine which markers to use as covariates
  
  if(n.perm > 0){
    if(verbose){cat("\n\nPerforming permutations to calculate significance threshold...\n")}			
    singlescan.obj$alpha.thresh <- genome.wide.threshold.1D(data.obj, geno.obj, n.perm = n.perm, scan.what = scan.what, ref.allele = ref.allele, alpha = alpha, model.family = model.family, n.cores = n.cores, run.parallel = TRUE)
  }else{
    cat("Not performing permutations...\n")	
  }
  
  
  #check for a reference allele, pull it out of the 
  #allele names here and add it to the data object
  if(length(ref.allele) != 1){ #add a check for the reference allele
    stop("You must specify one reference allele")
  }
  ref.col <- which(dimnames(gene)[[allele.dim]] == ref.allele)
  if(length(ref.col) == 0){
    stop("I can't find reference allele: ", ref.allele)
  }
  new.allele.names <- dimnames(gene)[[allele.dim]][-ref.col]
  singlescan.obj$ref.allele <- ref.allele	
  
  #=====================================================================
  #internal functions
  #=====================================================================
  
  #This function takes the results from get.stats.multiallele
  #and parses them into the final arrays
  add.results.to.array <- function(result.array, results.list, stat.name, is.covar = FALSE){
    row.num <- which(rownames(results.list[[1]][[1]]) == stat.name)
    #find the next spot to put data
    #Each successive phenotype is stored
    #in the 2nd dimension of the array
    next.spot.locale <- min(which(is.na(result.array[nrow(result.array),,1])))
    result.mat <- t(sapply(results.list, function(x) as.vector(x[[1]][row.num,])))
    
    if(is.covar){ #if we are looking at a covariate, we need to expand it
      result.mat <- matrix(result.mat, nrow = ncol(result.mat), ncol = dim(result.array)[3])
    }
    
    placement.locale <- match(names(results.list), rownames(result.array))
    # na.locale <- which(is.na(result.array[,next.spot.locale,1]))
    result.array[placement.locale,next.spot.locale,] <- result.mat
    return(result.array)	
  }
  
  add.flat.results.to.array <- function(result.array, model, pheno.num){
    betas <- as.vector(model[[2]]$beta)
    num.markers <- dim(result.array)[1]
    num.alleles <- dim(result.array)[3]
    start.pos <- 1
    for(i in 1:num.markers){
      locus.results <- betas[start.pos:(start.pos+num.alleles-1)]
      result.array[i,pheno.num,] <- locus.results
      start.pos = start.pos + num.alleles
    }
    return(result.array)
  }
  
  
  #=====================================================================
  #end of internal functions
  #=====================================================================
  
  
  #=====================================================================
  #begin code for multi-allelic cross
  #=====================================================================
  #In the multi-allele case, we want to collect
  #three 3D arrays each of num.marker by num.pheno by num.allele:
  #array of t statistics (for plotting p vals of regressions)
  #array of effects (betas) (for effect plots)
  #array of covar flags (for use in pair.scan)
  
  t.stat.array <- array(dim = c(dim(gene)[[locus.dim]]+n.covar, dim(pheno)[2], (dim(gene)[[allele.dim]]-1)))
  effect.array <- array(dim = c(dim(gene)[[locus.dim]]+n.covar, dim(pheno)[2], (dim(gene)[[allele.dim]]-1)))
  dimnames(t.stat.array) <- dimnames(effect.array) <- list(c(dimnames(gene)[[locus.dim]], covar.names), dimnames(pheno)[[2]], new.allele.names)
  
  #make arrays to hold locus-by-locus stats
  locus.score.scores <- matrix(NA, ncol = n.phe, nrow = dim(gene)[[locus.dim]]+n.covar)
  locus.p.values <- matrix(NA, ncol = n.phe, nrow = dim(gene)[[locus.dim]]+n.covar)
  rownames(locus.score.scores) <- rownames(locus.p.values) <- c(dimnames(gene)[[locus.dim]], covar.names)
  colnames(locus.score.scores) <- colnames(locus.p.values) <- colnames(pheno)
  
  for(i in 1:n.phe){
    if(verbose){cat("\nScanning trait:", colnames(pheno)[i], "\n")}
    #take out the response variable
    phenotype <- matrix(pheno[,i], ncol = 1)
    ph.family = model.family[i]
    
    #get corrected genotype and phenotype values for each phenotype-chromosome pair
    if(use.kinship){
      sink(file.path(data.obj$results_path,"regress.warnings")) #create a temporary output file for the regress warnings
      # TODO check if dim(kin.obj)[1] == length(phenoV) == length(covarV) when using covariates
      cor.data <- lapply(chr.which, function(x) kinship.on.the.fly(kin.obj, gene, chr1 = x, chr2 = x, phenoV = phenotype, covarV = covar.table))
      sink(NULL)
      
      names(cor.data) <- chr.which
    }else{
      cor.data <- vector(mode = "list", length = 1)
      cor.data[[1]] <- list("corrected.pheno" = phenotype, "corrected.geno" = gene, "corrected.covar" = covar.table)
    }
    
    results.by.chr <- vector(mode = "list", length = length(cor.data))							
    
    for(ch in 1:length(cor.data)){
      if(use.kinship){cat(" Chr", ch, "... ", sep = "")}
      if(length(cor.data) == 1){chr.locale <- 1:dim(gene)[3]}
      if(length(cor.data) > 1){chr.locale <- which(data.obj$chromosome == names(cor.data)[ch])}
      c.geno <- cor.data[[ch]]$corrected.geno[,,chr.locale,drop=FALSE]
      c.pheno <- cor.data[[ch]]$corrected.pheno
      c.covar <- cor.data[[ch]]$corrected.covar
      
      if (run.parallel) {
        
        cl <- parallel::makeCluster(n.cores)
        doParallel::registerDoParallel(cl)
        # the following line adds package variables to the parallel worker environments
        parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
        
        parallel::clusterEvalQ(cl, library(cape))
        
        # Copy functions in the workspace to the workers
        # funcs <- as.vector(utils::lsf.str(envir=.GlobalEnv))
        # parallel::clusterExport(cl, funcs, envir=.GlobalEnv)
        # cat("1\n")
        # parallel::clusterExport(cl, "check.geno", envir=.GlobalEnv)
        # cat("2\n")?
        results.by.chr <- foreach::foreach(x = 1:dim(c.geno)[locus.dim], .packages = 'cape') %dopar% {
          # Note that "Show Diagnostics" in RStudio will throw a warning that the `x` variable below is undefined
          # but it actually is defined in the foreach line above. You can safely ignore the warning.
          cape::get.stats.multiallele(phenotype = c.pheno, genotype = c.geno[,,x], covar.table = c.covar, ph.family, ref.col)
        }
        parallel::stopCluster(cl)
        
      } else {
        
        results.by.chr <- c()
        index <- 1:dim(c.geno)[locus.dim]
        for (x in index) {
          results.by.chr[[x]] <- get.stats.multiallele(phenotype = c.pheno, genotype = c.geno[,,x], covar.table = c.covar, ph.family, ref.col)
        }
        
      }
      names(results.by.chr) <- dimnames(c.geno)[[locus.dim]]	
      
      t.stat.array <- add.results.to.array(result.array = t.stat.array, results.list = results.by.chr, stat.name = "t.stat")
      effect.array <- add.results.to.array(effect.array, results.by.chr, "slope")
      locus.score.scores[chr.locale,i] <- unlist(lapply(results.by.chr, function(x) x$score))
      
    } #end looping through data corrected by chromosome (loco)
    
    #if there are covariates, run them through too
    if(!is.null(covar.table)){
      cat("\n")
      #get corrected genotype and phenotype values for the overall kinship matrix
      if(use.kinship){
        sink(file.path(data.obj$results_path,"regress.warnings"))
        # TODO check if dim(kin.obj)[1] == length(phenoV) == length(covarV) when using covariates
        cor.data <- kinship.on.the.fly(kin.obj, gene, phenoV = phenotype, covarV = covar.table)
        sink(NULL) #stop sinking to the file
      }else{
        cor.data <- list("corrected.pheno" = phenotype, "corrected.geno" = gene, "corrected.covar" = covar.table)
      }
      c.geno <- cor.data$corrected.geno
      c.pheno <- cor.data$corrected.pheno
      c.covar <- cor.data$corrected.covar 
      covar.results <- apply(c.covar, 2, function(x) get.stats.multiallele(c.pheno, x, c.covar, ph.family, ref.col))
      names(covar.results) <- data.obj$p_covar
      t.stat.array <- add.results.to.array(result.array = t.stat.array, results.list = covar.results, stat.name = "t.stat", is.covar = TRUE)
      effect.array <- add.results.to.array(effect.array, covar.results, "slope", is.covar = TRUE)
      first.na <- min(which(is.na(locus.score.scores[,i])))
      locus.score.scores[(first.na):(first.na+n.covar-1),i] <- unlist(lapply(covar.results, function(x) x$score))
    } #end case for an existing covariate table
    
  }      #end looping through phenotypes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  if(verbose){cat("\n")}
  
  singlescan.obj$singlescan.effects <- effect.array
  singlescan.obj$singlescan.t.stats <- t.stat.array
  # singlescan.obj$locus.p.vals <- locus.p.values
  singlescan.obj$locus.score.scores <- locus.score.scores
  singlescan.obj$n.perm <- n.perm
  
  unlink("regress.warnings") #remove the temporary file
  return(singlescan.obj)
  
}

