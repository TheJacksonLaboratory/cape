#' Select markers for pairscan randomly with uniform distribution
#'
#' This function selects markers for the pairscan
#' uniformly at random. I'm hoping that this minimizes
#' the skew in the m12/m21 null distribution. 
#' It has the option of forcing
#' specific markers to be included for replication of
#' previous results, etc.
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#' @param required.markers optional, a list of markers that must be included
#' @param num.alleles default = 500
#' @param verbose boolean, default = FALSE
#' 
select.markers.for.pairscan.uniform <- function(data.obj, geno.obj, required.markers = NULL, num.alleles = 500, verbose = FALSE){
  
  #require(abind)
  data.obj$marker_selection_method <- "uniform"
  ref.allele <- data.obj$ref_allele
  
  geno <- get.geno(data.obj, geno.obj)
  if(num.alleles >= dim(geno)[3]){
    required.markers = dimnames(geno)[[3]]
  }
  
  data.obj$ref.allele <- ref.allele
  
  other.allele <- setdiff(dimnames(geno)[[2]], ref.allele)
  
  if(length(required.markers) > 0){
    selected.idx <- sort(match(required.markers, dimnames(geno)[[3]]))
  }else{
    selected.idx <- NULL
  }
  
  #if the required markers are only some of the markers we want, 
  #select some more
  if(num.alleles > length(required.markers)){
    #select markers uniformly at random
    selected.idx <- sort(c(selected.idx, round(runif((num.alleles-length(required.markers)), 1, dim(geno)[[3]]))))
  }
  
  geno.for.pairscan <- geno[,other.allele,selected.idx]
  
  if(verbose){cat("Removing markers that are not linearly independent...\n")}
  data.obj$geno_for_pairscan <- geno.for.pairscan
  geno.ind <- get.linearly.independent(data.obj)
  if(verbose){
    cat(length(geno.ind[[2]]), "allele(s) rejected.\n")
    cat("Final alleles selected:", "\t", ncol(geno.ind$independent_markers), "\n")
  }
  
  rownames(geno.ind$independent.markers) <- rownames(data.obj$pheno)
  data.obj$geno_for_pairscan <- geno.ind$independent.markers		
  data.obj$marker.selection.method = "uniform"
  
  return(data.obj)
}