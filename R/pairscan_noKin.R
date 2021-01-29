#' Perform pairscan without a kinship correction
#' 
#' This internal function is called by \code{\link{pairscan}}
#' when no kinship correction is requested. It can be
#' compared to \code{\link{pairscan_kin}}. 
#' It fits pairwise linear models to estimate the effects of 
#' marker pairs on each trait. 
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param pheno_mat The matrix of trait values with individuals in rows
#' and traits in columns.
#' @param geno_mat The matrix of genotypes to be tested
#' @param covar_table The matrix of covariates with individuals in rows.
#' @param marker_pairs A two-column matrix containing the marker pairs
#' to be tested.
#' @param n_perm The number of permutations to be run
#' @param verbose A logical value indicating whether to 
#' print progress to the screen
#' @param run_parallel A logical value indicating 
#' whether parallel processing should be used
#' @param n_cores The number of cores to be used if run_parallel is TRUE
#' 
#' @return This function calls \code{\link{one_pairscan_parallel}} and
#' returns results for each trait as an element in a list.
#' @keywords internal
#' 
pairscan_noKin <- function(data_obj, pheno_mat, geno_mat, covar_table, 
                           marker_pairs, n_perm, verbose = FALSE, run_parallel = FALSE, 
                           n_cores = NULL){
  
  num_pheno <- dim(pheno_mat)[2]
  # cat("num_pheno:", num_pheno, "\n")
  
  results_list <- vector(mode = "list", length = num_pheno)
  names(results_list) <- colnames(pheno_mat)
  
  
  for(p in 1:num_pheno){ 
    if(verbose){
      cat("\nScanning phenotype ", colnames(pheno_mat)[p], ":\n", sep = "")
    }
    pairscan_results <- one_pairscan_parallel(data_obj, phenotype_vector = pheno_mat[,p], 
    		genotype_matrix = geno_mat, covar_vector = covar_table, paired_markers = marker_pairs, 
    		n_perm = 0, verbose = verbose, n_cores = n_cores, run_parallel = run_parallel)
    results_list[[p]] <- pairscan_results[[1]]
  } #end looping over phenotypes
  
  return(results_list)
  
  
}