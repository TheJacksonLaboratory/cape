#' Selects markers for pairscan based on effect sizes
#' 
#' This function takes in a singlescan object
#' and selects markers for the pairscan based on
#' the effect sizes of the markers originally selected
#' for the pairscan
#' Currently, this function is only used in generating a null
#' distribution, not in selecting markers initially. So it is
#' essentially not used anymore. At some point, I'd like to use
#' is for initial marker selection as well.
#' TODO !!!This function does not yet handle multi-parent crosses!!!
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param singlescan.obj a single scan object
#' @param geno.obj a genotype object
#' @param verbose boolean, default = FALSE
#' 
select.markers.for.pairscan.dist <- function(data.obj, singlescan.obj, geno.obj, verbose = FALSE){
  
  require(abind)
  
  if(is.null(data.obj$marker_selection_method)){
    data.obj$marker_selection_method <- "effects.dist"
  }
  
  ref.allele <- data.obj$ref_allele
  
  #===============================================================
  # get the distribution of effect sizes across the traits
  # we want to match this distribution
  #===============================================================
  geno.for.pairscan <- data.obj$geno.for.pairscan
  split.markers <- strsplit(colnames(geno.for.pairscan), "_")
  just.markers <- unlist(lapply(split.markers, function(x) x[1]))
  just.alleles <- unlist(lapply(split.markers, function(x) x[2]))
  
  if(class(singlescan.obj) == "list"){ 
    results <- abs(singlescan.obj$singlescan.t.stats) #an actual singlescan object
  }else{
    results <- abs(singlescan.obj) #a singlescan matrix for calculating pairscan null distribution
  }
  
  covar.info <- get.covar(data.obj)
  results.no.covar <- results[which(!rownames(results) %in% covar.info$covar.names),,,drop=FALSE]
  
  selected.effects <- t(mapply(function(x, y) results.no.covar[x,,y,drop=FALSE], just.markers, just.alleles))
  if(class(singlescan.obj) == "list"){
    colnames(selected.effects) <- dimnames(singlescan.obj$singlescan.effects)[[2]]
  }else{
    colnames(selected.effects) <- dimnames(singlescan.obj)[[2]]	
  }
  
  orig.effects <- selected.effects	
  
  
  #===============================================================
  # find all effects in the current single scan
  #===============================================================
  
  if(class(singlescan.obj) == "list"){
    ref.allele <- singlescan.obj$ref.allele
    data.obj$ref_allele <- ref.allele
  }else{
    ref.allele <- data.obj$ref_allele
  }
  
  selected.markers <- vector(mode = "list", length = ncol(orig.effects))
  # selected.effects <- vector(mode = "list", length = ncol(orig.effects))
  for(i in 1:length(selected.markers)){
    # TODO find the function def for get.nearest.pt()
    selected.markers[[i]] <- unlist(lapply(orig.effects[,i], function(x) rownames(results.no.covar)[get.nearest.pt(results.no.covar[,i,], x)]))
  }
  u_selected <- unique(unlist(selected.markers))		
  
  # plot(density(orig.effects))
  # points(density(unlist(selected.effects)), col = "red", type = "l")
  
  geno <- get.geno(data.obj, geno.obj)
  marker.locale <- sort(match(u_selected, dimnames(geno)[[3]]))
  
  geno.for.pairscan <- geno[,ref.allele,marker.locale]
  
  #===============================================================
  # check the final list for linear independence
  #===============================================================
  
  if(verbose){cat("Checking for linear independence...\n")}
  data.obj$geno.for.pairscan <- geno.for.pairscan
  geno.ind <- get.linearly.independent(data.obj)
  
  rownames(geno.ind$independent.markers) <- rownames(data.obj$pheno)
  data.obj$geno_for_pairscan <- geno.ind$independent.markers
  
  if(verbose){
    cat(length(geno.ind[[2]]), "allele(s) rejected.\n")
    cat("Final alleles selected:", "\t", ncol(geno.ind$independent.markers), "\n")
  }
  data.obj$marker_selection_method = "effects.dist"				
  return(data.obj)
  
  
}