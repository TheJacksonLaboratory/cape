#' This function returns the marker chromosome number with covariates 
#' returning Chr 0
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param markers covariate names
#' 
get.marker.chr <- function(data.obj, markers, character.names = TRUE){
  
  und.check <- grep("_", markers)
  if(length(und.check) > 0){
    markers <- sapply(strsplit(markers, "_"), function(x) x[1])
  }
  
  if(is.null(character.names)){
    is.char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
  }
  if(character.names){
    is.char <- TRUE
  }
  if(!character.names){
    is.char <- FALSE
  }
  
  
  if(is.char){
    marker.chr <- data.obj$chromosome[match(markers, data.obj$geno.names[[3]])]
    na.locale <- which(is.na(marker.chr))
    if(length(na.locale) > 0){
      covar.info <- get.covar(data.obj)
      geno.covar <- which(covar.info$covar.type == "g")
      if(length(geno.covar) > 0){
        unassigned <- markers[na.locale]
        unassigned.locale <- match(unassigned, data.obj$g.covar[1,])
        geno.covar.chr <- data.obj$g.covar[2,unassigned.locale]
        marker.chr[na.locale] <- geno.covar.chr
      }
    }
    #The rest are phenotypic covariates which are assigned chr 0
    na.locale <- which(is.na(marker.chr))
    if(length(na.locale) > 0){
      marker.chr[na.locale] <- 0
    }
    return(marker.chr)
  }
  
  
  
  if(!is.char){
    marker.chr <- data.obj$chromosome[match(markers, data.obj$marker.num)]
    na.locale <- which(is.na(marker.chr))
    if(length(na.locale) > 0){
      covar.info <- get.covar(data.obj)
      geno.covar <- which(covar.info$covar.type == "g")
      if(length(geno.covar) > 0){
        unassigned <- markers[na.locale]
        unassigned.locale <- match(unassigned, colnames(data.obj$g.covar))
        geno.covar.chr <- data.obj$g.covar[2,unassigned.locale]
        marker.chr[na.locale] <- geno.covar.chr
      }
    }
    #The rest are phenotypic covariates which are assigned chr 0
    na.locale <- which(is.na(marker.chr))
    if(length(na.locale) > 0){
      marker.chr[na.locale] <- 0
    }
    return(marker.chr)
  }
  
  
  
  
  
  
}
