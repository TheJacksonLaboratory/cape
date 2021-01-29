#' Removes markers from data_obj that are not present in the geno_obj
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#'
#' @return The data_obj is returned, and any markers that were
#' not present in geno_obj are removed from data_obj$geno_names
#'
#' 
#' @keywords internal
#' 
compare_markers <- function(data_obj, geno_obj){	
  geno <- get_geno(data_obj, geno_obj)
  missing_markers <- setdiff(data_obj$geno_names[[3]], dimnames(geno)[[3]])
  if(length(missing_markers) > 0){
    message("Removing markers from data_obj that are not present in the geno_obj\n")
    data_obj <- remove_markers(data_obj, missing_markers)
  }
  return(data_obj)	
}
