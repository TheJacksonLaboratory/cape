#' Runs the CAPE algorithm
#'
#' This function assumes you already have all required libraries and functions loaded.
#'
#' @param cape.obj the S4 class from \code{\link{Cape}}
#' @param parameter.file a full path string to the YAML file containing configuration parameters
#  TODO the config parameters together with the parameters for the cape.obj should be sufficient to reproduce a run of CAPE
#' @param results.dir a full path string to an existing empty directory. An error is thrown if the directory is not empty.
#' @param verbose boolean, output goes to stdout
#' @param run.parallel boolean, if TRUE runs certain parts of the code as parallel blocks
#'
#' @return None, output artifacts are saved to the results.dir directory
#'
#' @export
compare.markers <- function(data.obj, geno.obj){	
  geno <- get.geno(data.obj, geno.obj)
  missing.markers <- setdiff(data.obj$geno.names[[3]], dimnames(geno)[[3]])
  if(length(missing.markers) > 0){
    message("Removing markers from data.obj that are not present in the geno.obj")
    data.obj <- remove.markers(data.obj, missing.markers)
  }
  return(data.obj)	
}
