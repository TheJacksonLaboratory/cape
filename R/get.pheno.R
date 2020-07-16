#' Get the phenotype matrix
#' 
#' This function can return a number of different trait matrices
#' depending on the arguments.
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param scan.what A character string. One of "eigentraits", "normalized.trait",
#' or "raw.traits." If "eigentraits" the function returns the eigentraits matrix.
#' If "normalized.traits" the function returns the trait matrix after mean-centering
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
#' @export
get.pheno <- function(data.obj, scan.what = c("eigentraits", "normalized.traits", "raw.traits"), covar = NULL){
  
  #If the user does not specify a scan_what in the data.obj,
  #default to eigentraits, basically, if eigen,
  #et, or ET are anywhere in the string, use the
  #eigentraits, otherwise, use raw phenotypes
  if (is.null(data.obj$scan_what)) {
    data.obj$scan_what <- "eigentraits"
  }
  scan.what <- data.obj$scan_what
  is.ET <- c(grep("eig", scan.what, ignore.case = TRUE), grep("ET", scan.what, ignore.case = TRUE))
  is.raw <- grep("w", scan.what, ignore.case = TRUE)	 
  is.norm <- grep("o", scan.what, ignore.case = TRUE)	 
  
  
  if(length(is.ET) > 0){ 
    el.idx <- "ET"
    if(length(el.idx) == 0){
      stop("There are no eigentraits. Run get.eigentraits() to generate eigentraits.")
    }
  }
  
  if(length(is.raw) > 0){
    el.idx <- "pheno"
  }
  
  if(length(is.norm) > 0){
    el.idx <- "raw_pheno"
    if(length(el.idx) == 0){
      el.idx <- "pheno"
    }
  }
  
  pheno <- data.obj[[el.idx]]
  
  if(!is.null(covar)){
    models <- apply(pheno, 2, function(x) lm(x~covar.info$covar.table[,covar.locale,drop=FALSE]))
    resids <- lapply(models, residuals)
    resid.table <- t(matrix(unlist(resids, use.names = FALSE), byrow = TRUE))
    return(resid.table)
  }
  
  return(pheno)
  
  
}