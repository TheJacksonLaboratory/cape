#' Get marker names
#' 
#' Given a vector of marker numbers, this function 
#' returns the name of each marker.
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param markers A vector of marker numbers 
#' 
#' @return A vector the same length as the input 
#' markers vector indicating the name of each marker
#' 
#' @examples 
#' \dontrun{
#' marker_names <- get_marker_chr(data_obj, 1:10)
#' }
#'
#' @export
get_marker_name <- function(data_obj, markers){
  
  if(is.null(markers)){return(NULL)}
  
  und_check <- grep("_", markers[1])
  if(length(und_check) > 0){
    markers <- sapply(strsplit(markers, "_"), function(x) x[1])
  }
  
  is_char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
  marker_names <- data_obj$geno_names[[3]]
  
  
  if(is_char){
    #use the marker_names vector first to translate
    marker_name <- marker_names[match(markers, marker_names)]
    
    #if there are any markers we didn't translate, look in the 
    #covariate tables for marker numbers
    na_locale <- which(is.na(marker_name))
    
    if(length(na_locale) > 0){
      covar_info <- get_covar(data_obj)
      marker_name[na_locale] <- covar_info$covar_names[match(markers[na_locale], covar_info$covar_names)]
    }
  }else{
    #use the marker_names vector first to translate
    marker_name <- marker_names[match(markers, data_obj$marker_num)]
    
    #if there are any markers we didn't translate, look in the 
    #covariate tables for marker numbers
    na_locale <- which(is.na(marker_name))
    
    if(length(na_locale) > 0){
      covar_info <- get_covar(data_obj)
      marker_name[na_locale] <- covar_info$covar_names[match(markers[na_locale], colnames(covar_info$covar_table))]
    }
  }
  
  return(marker_name)
}