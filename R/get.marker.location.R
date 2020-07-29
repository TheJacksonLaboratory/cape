#' Get marker genomic position
#' 
#' Given a vector of marker names or numers, this 
#' function returns the genomic coordinates for
#' each marker, not including the chromosome number,
#' which is retrieved using \link{\code{get.marker.chr}}.
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param markers A vector of marker names 
#' 
#' @return A vector the same length as the input markers vector
#' indicating the genomic coordinate of each marker.
#'
#' @export

get.marker.location <- function(data.obj, markers){
  
  
  markers.with.na <- markers
  not.na.pos <- which(!is.na(markers.with.na))
  markers <- markers[not.na.pos]
  
  
  und.check <- grep("_", markers)
  if(length(und.check) > 0){
    markers <- sapply(strsplit(markers, "_"), function(x) x[1])
  }
  
  is.char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
  
  if(is.char){
    marker.loc <- data.obj$marker_location[match(markers, data.obj$geno_names[[3]])]
    na.locale <- which(is.na(marker.loc))
    if(length(na.locale) > 0){
      covar.info <- get.covar(data.obj)
      geno.covar <- which(covar.info$covar.type == "g")
      if(length(geno.covar) > 0){
        geno.covar.locale <- match(markers[na.locale], covar.info$covar.names[geno.covar])
        geno.covar.loc <- data.obj$g_covar[3,geno.covar.locale]
        marker.loc[na.locale] <- geno.covar.loc
      }
    }
    na.locale <- which(is.na(marker.loc))
    #if there are still missing values, these are 
    #phenotypic covariates. They get dummy positions
    if(length(na.locale) > 0){ 
      pheno.covar <- which(covar.info$covar.type == "p")
      if(length(pheno.covar) > 0){
        marker.loc[na.locale] <- 1:length(na.locale)
      }
    } #end case for if there are still missing markers
    final.marker.pos <- rep(NA, length(markers.with.na))
    final.marker.pos[not.na.pos] <- marker.loc
    return(final.marker.pos)
  }
  
  
  if(!is.char){
    marker.loc <- data.obj$marker_location[match(markers, data.obj$marker_num)]
    na.locale <- which(is.na(marker.loc))
    if(length(na.locale) > 0){
      covar.info <- get.covar(data.obj)
      covar.locale <- which(!is.na(match(colnames(covar.info$covar.table), markers)))
      if(length(covar.locale) == 0){
        marker.loc[na.locale] <- match(markers[na.locale], covar.info$covar.names)
      }else{
        marker.loc[na.locale] <- match(markers[na.locale], colnames(covar.info$covar.table))
      }
    }
    final.marker.pos <- rep(NA, length(markers.with.na))
    final.marker.pos[not.na.pos] <- marker.loc
    return(final.marker.pos)
  }
  
}
