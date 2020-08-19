#' Checks for underscores in marker names
#' 
#' This function looks for underscores in 
#' marker names. Beause cape uses underscores
#' to separate marker from allele names, we
#' cannot have any in the marker names
#'
#' @param data.obj a \code{\link{Cape}} object

check.underscore <- function(data.obj){
  
  marker.names <- data.obj$geno_names[[3]]
  under.locale <- grep("_", marker.names)
  if(length(under.locale) > 0){
    stop("Underscores have been detected in some marker names.\nPlease use delete_underscore() to remove these before proceeding.")
  }
}

