#' Perform pairscan without a kinship correction
#' 
#' This internal function is called by \link{\code{pairscan}}
#' when no kinship correction is requested. It can be
#' compared to \link{\code{pairscan.kin}}. 
#' It fits pairwise linear models to estimate the effects of 
#' marker pairs on each trait. 
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param pheno.mat The matrix of trait values with individuals in rows
#' and traits in columns.
#' @param geno.mat The matrix of genotypes to be tested
#' @param covar.table The matrix of covariates with individuals in rows.
#' @param marker.pairs A two-column matrix containing the marker pairs
#' to be tested.
#' n.perm The number of permutations to be run
#' @param verbose A logical value indicating whether to 
#' print progress to the screen
#' @param run.parallel A logical value indicating 
#' whether parallel processing should be used
#' @param n.cores The number of cores to be used if run.parallel is TRUE
#' 
#' @return This function calls \link{\code{one.pairscan.parallel}} and
#' returns results for each trait as an element in a list.
#' 
pairscan.noKin <- function(data.obj, pheno.mat, geno.mat, covar.table, marker.pairs, n.perm, verbose = FALSE, run.parallel = FALSE, n.cores = NULL){
  
  num.pheno <- dim(pheno.mat)[2]
  # cat("num.pheno:", num.pheno, "\n")
  
  results.list <- vector(mode = "list", length = num.pheno)
  names(results.list) <- colnames(pheno.mat)
  
  
  for(p in 1:num.pheno){ 
    if(verbose){
      cat("\nScanning phenotype ", colnames(pheno.mat)[p], ":\n", sep = "")
    }
    pairscan.results <- one.pairscan.parallel(data.obj, phenotype.vector = pheno.mat[,p], 
    		genotype.matrix = geno.mat, covar.vector = covar.table, paired.markers = marker.pairs, 
    		n.perm = 0, verbose = verbose, n.cores = n.cores, run.parallel = run.parallel)
    results.list[[p]] <- pairscan.results[[1]]
  } #end looping over phenotypes
  
  return(results.list)
  
  
}