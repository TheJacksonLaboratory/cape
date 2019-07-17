#' This function returns the marker names given the numbers
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param markers covariate names
#' 
get.marker.name <- function(data.obj, markers){
  
  if(is.null(markers)){return(NULL)}
  
  und.check <- grep("_", markers[1])
  if(length(und.check) > 0){
    markers <- sapply(strsplit(markers, "_"), function(x) x[1])
  }
  
  is.char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
  marker.names <- data.obj$geno_names[[3]]
  
  
  if(is.char){
    #use the marker.names vector first to translate
    marker.name <- marker.names[match(markers, marker.names)]
    
    #if there are any markers we didn't translate, look in the 
    #covariate tables for marker numbers
    na.locale <- which(is.na(marker.name))
    
    if(length(na.locale) > 0){
      covar.info <- get.covar(data.obj)
      marker.name[na.locale] <- covar.info$covar.names[match(markers[na.locale], covar.info$covar.names)]
    }
  }else{
    #use the marker.names vector first to translate
    marker.name <- marker.names[match(markers, data.obj$marker.num)]
    
    #if there are any markers we didn't translate, look in the 
    #covariate tables for marker numbers
    na.locale <- which(is.na(marker.name))
    
    if(length(na.locale) > 0){
      covar.info <- get.covar(data.obj)
      marker.name[na.locale] <- covar.info$covar.names[match(markers[na.locale], colnames(covar.info$covar.table))]
    }
  }
  
  return(marker.name)
  
  
  
}