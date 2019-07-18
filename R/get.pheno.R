#' Finds alleles from each founder were seleted for the pairscan
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param color.scheme string "CC/DO" or "other"
#' @param pdf.label string
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