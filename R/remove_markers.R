#' Removes genetic markers
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param markers_to_remove A vector of marker names to be removed.
#' 
#' @examples 
#' \dontrun{
#' #remove markers on chromosome 1
#' marker_idx <- which(data_obj$chromosome == 1)
#' data_obj <- remove_markers(data_obj, marker_idx)
#' }
#' 
#' @export
remove_markers <- function(data_obj, markers_to_remove){
  
  marker_idx <- match(markers_to_remove, data_obj$geno_names[[3]])
  if(any(is.na(marker_idx))){marker_idx <- markers_to_remove}
  
  #take out the markers from the data_obj meta data	
  data_obj$chromosome <- data_obj$chromosome[-marker_idx]
  data_obj$marker_num <- data_obj$marker_num[-marker_idx]
  data_obj$marker_location <- data_obj$marker_location[-marker_idx]
  data_obj$geno_names[[3]] <- data_obj$geno_names[[3]][-marker_idx]
  
  return(data_obj)
}

