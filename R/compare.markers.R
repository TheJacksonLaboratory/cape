#' Removes markers from data.obj that are not present in the geno.obj
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object. If this is not supplied then it is generated here.
#'
#' @return \code{list("data.obj" = data.obj, "geno.obj" = geno.obj)}
#'
compare.markers <- function(data.obj, geno.obj){	
  geno <- get.geno(data.obj, geno.obj)
  missing.markers <- setdiff(data.obj$geno_names[[3]], dimnames(geno)[[3]])
  if(length(missing.markers) > 0){
    cat("Removing markers from data.obj that are not present in the geno.obj\n")
    data.obj <- remove.markers(data.obj, missing.markers)
  }
  return(data.obj)	
}
