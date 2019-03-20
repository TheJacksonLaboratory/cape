#' Checks for underscores in marker names
#' 
#' Internal function.
#'
#' @param data.obj a \code{\link{Cape}} object
check.underscore <- function(data.obj){
  
  marker.names <- data.obj$geno_names[[3]]
  under.locale <- grep("_", marker.names)
  if(length(under.locale) > 0){
    stop("Underscores have been detected in some marker names.\nPlease use delete.underscore() to remove these before proceeding.")
  }
}

