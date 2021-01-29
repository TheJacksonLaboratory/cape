#' Get original indices for markers
#' 
#' Given a vector of marker names or numbers, this 
#' function returns the indices in the genotype 
#' matrix for each marker.
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param markers A vector of marker names 
#' 
#' @return A vector the same length as the input markers vector
#' indicating the index of each marker
#' @keywords internal
#'

get_marker_idx <- function(data_obj, markers){
  
  und_check <- grep("_", markers[1])
  if(length(und_check) > 0){
    markers <- sapply(strsplit(markers, "_"), function(x) x[1])
  }
  
  is_char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
  
  if(is_char){
    marker_loc <- match(markers, data_obj$geno_names[[3]])
    return(marker_loc)
  }
  
  if(!is_char){
    marker_loc <- match(markers, data_obj$marker_num)
    return(marker_loc)
  }
}
