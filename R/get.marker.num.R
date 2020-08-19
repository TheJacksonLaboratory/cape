#' Get numbers for markers
#' 
#' Given a vector of marker names or numbers, this 
#' function returns the numeric index of each marker.
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param markers A vector of marker names 
#' 
#' @return A vector the same length as the input markers vector
#' indicating the number of each chromosome.
#'

get.marker.num <- function(data.obj, markers){
  
  if(is.null(markers)){
    return(NULL)
  }
  
  und.check <- grep("_", markers[1])
  if(length(und.check) > 0){
    markers <- sapply(strsplit(markers, "_"), function(x) x[1])
  }
  
  is.char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
  
  
  if(is.char){
    #use the markers vector first to translate
    marker.num <- data.obj$marker_num[match(markers, data.obj$geno_names[[3]])]
    
    #if there are any markers we didn't translate, look in the 
    #covariate tables for marker numbers
    na.locale <- which(is.na(marker.num))
    
    if(length(na.locale) > 0){
      covar.info <- get.covar(data.obj)
      marker.num[na.locale] <- colnames(covar.info$covar.table)[match(markers[na.locale], covar.info$covar.names)]
    }
  }else{
    #use the markers vector first to translate
    marker.num <- data.obj$marker_num[match(markers, data.obj$marker_num)]
    
    #if there are any markers we didn't translate, look in the 
    #covariate tables for marker numbers
    na.locale <- which(is.na(marker.num))
    
    if(length(na.locale) > 0){
      covar.info <- get.covar(data.obj)
      marker.num[na.locale] <- colnames(covar.info$covar.table)[match(markers[na.locale], colnames(covar.info$covar.table))]
    }
  }
  
  
  return(marker.num)
}