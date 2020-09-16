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
get.pheno <- function(data.obj, 
scan.what = c("eigentraits", "normalized.traits", "raw.traits"), covar = NULL){
  

	scan.what <- scan.what[1] #default to eigentraits

  is.ET <- as.logical(length(c(grep("eig", scan.what, ignore.case = TRUE), grep("ET", scan.what, ignore.case = TRUE))))
  is.raw <- as.logical(length(grep("w", scan.what, ignore.case = TRUE)))
  is.norm <- as.logical(length(grep("o", scan.what, ignore.case = TRUE)))

  
  if(is.ET){ 
    el.idx <- "ET"
    if(is.null(data.obj$ET)){
      stop("There are no eigentraits. Run get.eigentraits() to generate eigentraits.")
    }
  }
  
  if(is.raw){
    el.idx <- "raw_pheno"
    if(is.null(data.obj$raw_pheno)){
    	el.idx <- "pheno"
    }
  }
  
  if(is.norm){
    el.idx <- "pheno"
  }
  
  pheno <- data.obj[[el.idx]]
  
  if(!is.null(covar)){
  	covar.info <- get.covar(data.obj)
  	covar.locale <- match(covar, covar.info$covar.names)
    models <- apply(pheno, 2, function(x) lm(x~covar.info$covar.table[,covar.locale,drop=FALSE]))
    resids <- lapply(models, residuals)
    resid.table <- Reduce("cbind", resids)
    colnames(resid.table) <- names(resids)
    return(resid.table)
  }
  
  return(pheno)
  
  
}