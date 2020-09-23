#' Mean-center and normalize phenotypes
#'
#' This function is a wrapper for mean-centering
#' normalizing and standardizing the trait matrix.
#' in a data_obj.
#'
#' @param data_obj a \code{\link{Cape}} object
#' mean_center a logical value indicating whether the
#' traits should be mean centered. If FALSE, the traits
#' are only normalized.
#' @param mean_center mean center
#'
#' @return the data object is returned. The pheno slot of
#' the data object will have normalized and/or mean-centered
#' traits. The function also preserves the original trait matrix
#' in a slot called raw_pheno.
#'
#'@export
norm_pheno <- function(data_obj, mean_center = TRUE){
  
  pheno <- data_obj$pheno
  raw_pheno <- data_obj$pheno #retain the raw phenotype values
  
  pheno <- apply(pheno, 2, rz_transform)
  
  if(mean_center){
    pheno <- apply(pheno, 2, center_std) #mean center and standardize the phenotypes
  }
  
  data_obj$pheno <- pheno #replace the raw phenotypes with scaled, normalized phenotypes (if we have done those things)
  data_obj$raw_pheno <- raw_pheno
  
  return(data_obj)
}