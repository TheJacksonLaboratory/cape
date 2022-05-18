#' Performs marker regression
#' 
#' This is an internal function called by \code{\link{pairscan_null}}
#' when generating the null distribution for significance testing. 
#' To perform permutations, we permute trait values, and then re-do
#' the singlescan, marker selection, and the pair scan on the permuted
#' data. This function performs the singlescan on the permuted data.
#' 
#' @param phenotype_vector A vector of phenotype values, one entry for each individual.
#' @param genotype_mat A matrix of genotype values with individuals in rows and markers 
#' in columns. Matrix entries contain the probability of the reference allele at each 
#' position for each individual.
#' @param model_family Indicates the model family of the phenotypes. This can be 
#'   either "gaussian" or "binomial".
#' @param ref_allele the reference allele from the Cape data_obj
#' @param covar_table A matrix of covariates with one row per individual.
#' @param run_parallel A logical value indicating whether multiple 
#' processors should be used
#' @param n_cores The number of processors to use if run_parallel is TRUE.
#' 
#' @return This function returns the t_statistics for all linear models 
#' testing the effects of each marker on the phenotype.
#' 
#' @import parallel
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @keywords internal
#'  
one_singlescanDO <- function(phenotype_vector, genotype_mat, model_family, ref_allele = "A", 
covar_table = NULL, run_parallel = FALSE, n_cores = 4){
  # declaring variable to prevent warning
  m <- NULL
  
  if(!run_parallel){n_cores = 1}
  
  gene <- genotype_mat
  
  #Get the dimension names to minimize confusion	
  geno_dims <- get_geno_dim()
  mouse_dim <- geno_dims[which(names(geno_dims) == "mouse")]
  allele_dim <- geno_dims[which(names(geno_dims) == "allele")]
  locus_dim <- geno_dims[which(names(geno_dims) == "locus")]
  
  
  ref_col <- which(dimnames(gene)[[allele_dim]] == ref_allele)
  new_allele_names <- dimnames(gene)[[allele_dim]][-ref_col]
  
  #=====================================================================
  #begin code for multi-allelic cross
  #=====================================================================
  
  #apply the modeling function to each marker column
  # !diagnostics suppress=m
  if (run_parallel) {
    
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    cape_dir_full <- find.package("cape")
    cape_dir <- gsub("cape_pkg/cape","cape_pkg", cape_dir_full)
    clusterExport(cl, "cape_dir", envir=environment())
    clusterEvalQ(cl, .libPaths(cape_dir))
    results <- foreach(m = 1:dim(gene)[[locus_dim]], .packages = 'cape', .export = c("get_stats_multiallele", "check_geno")) %dopar% {
      get_stats_multiallele(phenotype_vector, gene[,,m], covar_table = covar_table, 
      model_family, ref_col)
    }
    stopCluster(cl)
    
  } else {
    
    index <- 1:dim(gene)[[locus_dim]]
    results <- lapply(index, function(x) get_stats_multiallele(phenotype_vector, gene[,,x], 
      covar_table, model_family, ref_col))    
  }
  
  t_stat_array <- matrix(unlist(lapply(results, function(x) x[[1]]["t_stat",])), 
  ncol = length(new_allele_names), byrow = TRUE)
  colnames(t_stat_array) <- new_allele_names
  
  
  return(t_stat_array)		
}

