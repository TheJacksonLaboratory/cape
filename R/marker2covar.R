#' Allows users to change covariate assignments
#' 
#' This script allows users to change covariate assignments
#' manually.
#' markers <- rownames(covar.flags)[5:10]
#' pheno <- "ET1"
#' It adds the information to the data object
#' only the output from singlescan and pairscan
#' make separate objects.
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#' @param singlescan.obj a singlescan object is required if setting covariates by a t threshold.
#' @param covar.thresh 
#' @param markers
#'
#' @export
marker2covar <- function(data.obj, geno.obj = NULL, singlescan.obj = NULL, covar.thresh = NULL, markers = NULL){
  
  if(!is.null(covar.thresh)){
    oneD <- singlescan.obj$singlescan.results
    if(is.null(oneD)){stop("singlescan.obj is required if setting covariates by a t threshold.\n")}
  }
  
  
  geno.mat <- get.geno(data.obj, geno.obj)
  
  #if the user has specified a t threshold for covariate specification
  if(!is.null(covar.thresh)){
    
    marker.names <- data.obj$geno_names[[3]]
    
    singlescan.obj$covar_thresh <- covar.thresh
    
    covar.which <- lapply(oneD, function(x) which(x[,"t.stat"] >= covar.thresh))
    covar.names <- unique(unlist(lapply(covar.which, function(x) names(x))))
    browser()
    new.covar.locale <- get.marker.idx(data.obj, covar.names)
    new.covar <- geno.mat[,new.covar.locale,drop=FALSE]
    dimnames(new.covar)[[3]] <- covar.names
    
    snp.names <- get.marker.name(data.obj, covar.names)
    
    g.covar.info <- rbind(snp.names, data.obj$chromosome[new.covar.locale], data.obj$marker_location[new.covar.locale])
    colnames(g.covar.info) <- data.obj$marker_num[new.covar.locale]
    rownames(g.covar.info) <- c("name", "chromosome", "position")
    
    data.obj <- remove.markers(data.obj, markers.which = snp.names)
    data.obj$g_covar_table <- new.covar
    data.obj$g_covar <- g.covar.info
    return(data.obj)		
  } #end case for setting covariates by a threshold
  
  
  if(!is.null(markers)){
    marker.names <- data.obj$geno_names[[3]]
    
    marker.locale <- get.marker.idx(data.obj, markers)		
    new.covar <- geno.mat[,,marker.locale,drop=FALSE]
    ref.locale <- which(dimnames(new.covar)[[2]] == data.obj$ref_allele)
    the.rest <- setdiff(1:dim(new.covar)[[2]], ref.locale)
    covar.mat <- new.covar[,the.rest,]
    data.obj$g_covar_table <- covar.mat
    
    g.covar.info <- rbind(marker.names[marker.locale], data.obj$chromosome[marker.locale], data.obj$marker_location[marker.locale])
    colnames(g.covar.info) <- data.obj$marker_num[marker.locale]
    rownames(g.covar.info) <- c("name", "chromosome", "position")
    
    data.obj$g_covar <- g.covar.info
    data.obj <- remove.markers(data.obj, markers)
    
  }
  
  return(data.obj)
  
}