#' Removes unwanted markers
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param markers.to.remove
#' 
remove.markers <- function(data.obj, markers.to.remove){
  
  marker.idx <- match(markers.to.remove, data.obj$geno_names[[3]])
  if(any(is.na(marker.idx))){marker.idx <- markers.to.remove}
  
  #take out the markers from the data.obj meta data	
  data.obj$chromosome <- data.obj$chromosome[-marker.idx]
  data.obj$marker_num <- data.obj$marker_num[-marker.idx]
  data.obj$marker_location <- data.obj$marker_location[-marker.idx]
  data.obj$geno_names[[3]] <- data.obj$geno_names[[3]][-marker.idx]
  
  return(data.obj)
}

