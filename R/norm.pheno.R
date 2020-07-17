#' Mean-center and normalize phenotypes
#'
#' This function is a wrapper for mean-centering
#' normalizing and standardizing the trait matrix.
#' in a data.obj.
#'
#' @param data.obj a \link{\code{Cape}} object
#' mean.center a logical value indicating whether the
#' traits should be mean centered. If FALSE, the traits
#' are only normalized.
#'
#' @return the data object is returned. The pheno slot of
#' the data object will have normalized and/or mean-centered
#' traits. The function also preserves the original trait matrix
#' in a slot called raw.pheno.
#'

norm.pheno <- function(data.obj, mean.center = TRUE){
  
  pheno <- data.obj$pheno
  raw.pheno <- data.obj$pheno #retain the raw phenotype values
  
  pheno <- apply(pheno, 2, rz.transform)
  
  if(mean.center){
    pheno <- apply(pheno, 2, center.std) #mean center and standardize the phenotypes
  }
  
  data.obj$pheno <- pheno #replace the raw phenotypes with scaled, normalized phenotypes (if we have done those things)
  data.obj$raw_pheno <- raw.pheno
  
  return(data.obj)
}