#' Checks for underscores in marker names
#' 
#' This function looks for underscores in 
#' marker names. Beause cape uses underscores
#' to separate marker from allele names, we
#' cannot have any in the marker names
#'
#' @param data_obj a \code{\link{Cape}} object
#' @keywords internal

check_underscore <- function(data_obj){
  
  marker_names <- data_obj$geno_names[[3]]
  under_locale <- grep("_", marker_names)
  if(length(under_locale) > 0){
    stop("Underscores have been detected in some marker names.\nPlease use delete_underscore() to remove these before proceeding.")
  }
}

