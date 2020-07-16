#' Get original indices for markers
#' 
#' Given a vector of marker names or numbers, this 
#' function returns the indices in the genotype 
#' matrix for each marker.
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param markers A vector of marker names 
#' 
#' @return A vector the same length as the input markers vector
#' indicating the index of each marker
#'

get.marker.idx <- function(data.obj, markers){
  
  und.check <- grep("_", markers[1])
  if(length(und.check) > 0){
    markers <- sapply(strsplit(markers, "_"), function(x) x[1])
  }
  
  is.char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
  
  if(is.char){
    marker.loc <- match(markers, data.obj$geno_names[[3]])
    return(marker.loc)
  }
  
  if(!is.char){
    marker.loc <- match(markers, data.obj$marker_num)
    return(marker.loc)
  }
  
  
}
