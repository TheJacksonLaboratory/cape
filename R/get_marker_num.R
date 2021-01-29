#' Get numbers for markers
#' 
#' Given a vector of marker names or numbers, this 
#' function returns the numeric index of each marker.
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param markers A vector of marker names 
#' 
#' @return A vector the same length as the input markers vector
#' indicating the number of each chromosome.
#' @keywords internal
#'

get_marker_num <- function(data_obj, markers){
  
  if(is.null(markers)){
    return(NULL)
  }
  
  und_check <- grep("_", markers[1])
  if(length(und_check) > 0){
    markers <- sapply(strsplit(markers, "_"), function(x) x[1])
  }
  
  is_char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
  
  
  if(is_char){
    #use the markers vector first to translate
    marker_num <- data_obj$marker_num[match(markers, data_obj$geno_names[[3]])]
    
    #if there are any markers we didn't translate, look in the 
    #covariate tables for marker numbers
    na_locale <- which(is.na(marker_num))
    
    if(length(na_locale) > 0){
      covar_info <- get_covar(data_obj)
      marker_num[na_locale] <- colnames(covar_info$covar_table)[match(markers[na_locale], covar_info$covar_names)]
    }
  }else{
    #use the markers vector first to translate
    marker_num <- data_obj$marker_num[match(markers, data_obj$marker_num)]
    
    #if there are any markers we didn't translate, look in the 
    #covariate tables for marker numbers
    na_locale <- which(is.na(marker_num))
    
    if(length(na_locale) > 0){
      covar_info <- get_covar(data_obj)
      marker_num[na_locale] <- colnames(covar_info$covar_table)[match(markers[na_locale], colnames(covar_info$covar_table))]
    }
  }
  
  
  return(marker_num)
}