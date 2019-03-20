#This function returns the marker index for marker names or numbers

#' This function returns the marker index for marker names or numbers
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param markers
#'
#' @export
get.marker.idx <- function(data.obj, markers){
  
  und.check <- grep("_", markers[1])
  if(length(und.check) > 0){
    markers <- sapply(strsplit(markers, "_"), function(x) x[1])
  }
  
  is.char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
  
  if(is.char){
    marker.loc <- match(markers, data.obj$geno.names[[3]])
    return(marker.loc)
  }
  
  if(!is.char){
    marker.loc <- match(markers, data.obj$marker.num)
    return(marker.loc)
  }
  
  
}
