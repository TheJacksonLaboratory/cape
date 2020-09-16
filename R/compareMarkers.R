#' Removes markers from data.obj that are not present in the geno.obj
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#'
#' @return The data.obj is returned, and any markers that were
#' not present in geno.obj are removed from data.obj$geno_names
#'
#' @export
compareMarkers <- function(data.obj, geno.obj){	
  geno <- get.geno(data.obj, geno.obj)
  missing.markers <- setdiff(data.obj$geno_names[[3]], dimnames(geno)[[3]])
  if(length(missing.markers) > 0){
    cat("Removing markers from data.obj that are not present in the geno.obj\n")
    data.obj <- remove.markers(data.obj, missing.markers)
  }
  return(data.obj)	
}
