#' Get the phenotype matrix
#' 
#' This function can return a number of different trait matrices
#' depending on the arguments.
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param scan_what A character string. One of "eigentraits", "normalized.trait",
#' or "raw_traits." If "eigentraits" the function returns the eigentraits matrix.
#' If "normalized_traits" the function returns the trait matrix after mean-centering
#' and normalizing. If "raw.trait" the function returns the trait matrix before 
#' mean-centering and normalization were applied.
#' @param covar A character value indicating which, if any, covariates the traits
#' should be adjusted for. If covariates are specified, the function fits a linear
#' model to specify the traits with the covariates and returns the matrix of residuals
#' (i.e. the traits after adjusting for the covariates).
#'
#' @return A matrix in which each column is a trait, and each row is an individual.
#' The values correspond to the argument settings described above.
#' 
#' @importFrom stats residuals
#' 
#' @export
get_pheno <- function(data_obj, scan_what = c("eigentraits", "normalized_traits", "raw_traits"), covar = NULL){
  

	scan_what <- scan_what[1] #default to eigentraits

  is_ET <- as.logical(length(c(grep("eig", scan_what, ignore.case = TRUE), grep("ET", scan_what, ignore.case = TRUE))))
  is_raw <- as.logical(length(grep("w", scan_what, ignore.case = TRUE)))
  is_norm <- as.logical(length(grep("o", scan_what, ignore.case = TRUE)))

  
  if(is_ET){ 
    el_idx <- "ET"
    if(is.null(data_obj$ET)){
      stop("There are no eigentraits. Run get_eigentraits() to generate eigentraits.")
    }
  }
  
  if(is_raw){
    el_idx <- "raw_pheno"
    if(is.null(data_obj$raw_pheno)){
    	el_idx <- "pheno"
    }
  }
  
  if(is_norm){
    el_idx <- "pheno"
  }
  
  pheno <- data_obj[[el_idx]]
  
  if(!is.null(covar)){
  	covar_info <- get_covar(data_obj)
  	covar_locale <- match(covar, covar_info$covar_names)
    models <- apply(pheno, 2, function(x) lm(x~covar_info$covar_table[,covar_locale,drop=FALSE]))
    resids <- lapply(models, residuals)
    resid_table <- Reduce("cbind", resids)
    colnames(resid_table) <- names(resids)
    return(resid_table)
  }
  
  return(pheno)
}