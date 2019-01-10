#' Removes markers from cape.obj that are not present in the geno.obj
#'
#' @param cape.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object. If this is not supplied then it is generated here.
#'
#' @return \code{list("cape.obj" = cape.obj, "geno.obj" = geno.obj)}
#'
#' @export
compare.markers <- function(cape.obj, geno.obj){	
  geno <- get.geno(cape.obj, geno.obj)
  missing.markers <- setdiff(cape.obj$geno.names[[3]], dimnames(geno)[[3]])
  if(length(missing.markers) > 0){
    message("Removing markers from cape.obj that are not present in the geno.obj")
    cape.obj <- remove.markers(cape.obj, missing.markers)
  }
  return(cape.obj)	
}
