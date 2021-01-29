#' Delete underscores from marker names
#' 
#' Because cape uses underscores to separate marker names
#' from alleles, we need to delete any underscores that 
#' were present in the original marker names. This function
#' does that. 
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' 
#' @return lists containing \code{data_obj} and \code{geno_obj}
#' with underscores removed from \code{data_obj$geno_names}
#' and \code{dimnames(geno_obj)}. These two elements must be
#' entered separately into \code{\link{run_cape}}.
#' 
#' @keywords internal

delete_underscore <- function(data_obj, geno_obj = NULL){
  
  geno <- get_geno(data_obj, geno_obj)
  
  marker_names <- data_obj$geno_names[[3]]
  under_locale <- grep("_", marker_names)
  
  if(length(under_locale) > 0){
    bad_names <- marker_names[under_locale]
    new_names <- unlist(lapply(strsplit(bad_names, "_"), function(x) paste(x[1:length(x)], collapse = "")))
    
    data_obj$geno_names[[3]][under_locale] <- new_names
    dimnames(geno)[[3]][under_locale] <- new_names
    message("Removing underscores from marker names\n")
  }	
  
  results <- list(data_obj, geno)
  names(results) <- c("data_obj", "geno_obj")
  
  return(results)
  
}

