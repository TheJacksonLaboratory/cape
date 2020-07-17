#' Performs marker regression
#' 
#' This is an internal function called by \link{\code{pairscan.null}}
#' when generating the null distribution for significance testing. 
#' To perform permutations, we permute trait values, and then re-do
#' the singlescan, marker selection, and the pair scan on the permuted
#' data. This function performs the singlescan on the permuted data.
#' 
#' @param phenotype.vector A vector of phenotype values, one entry for each individual.
#' @param genotype.mat A matrix of genotype values with individuals in rows and markers 
#' in columns. Matrix entries contain the probability of the reference allele at each 
#' position for each individual.
#' @param model.family Indicates the model family of the phenotypes. This can be 
#'   either "gaussian" or "binomial".
#' @param ref.allele the reference allele from the Cape data.obj
#' @param covar.table A matrix of covariates with one row per individual.
#' @param run.parallel A logical value indicating whether multiple 
#' processors should be used
#' @param n.cores The number of processors to use if run.parallel is TRUE.
#' 
#' @return This function returns the t.statistics for all linear models 
#' testing the effects of each marker on the phenotype.


one.singlescanDO <- function(phenotype.vector, genotype.mat, model.family, ref.allele = "A", 
covar.table = NULL, run.parallel = FALSE, n.cores = 4){
  
  if(!run.parallel){n.cores = 1}
  
  gene <- genotype.mat
  
  #Get the dimension names to minimize confusion	
  geno.dims <- get_geno_dim()
  mouse.dim <- geno.dims[which(names(geno.dims) == "mouse")]
  allele.dim <- geno.dims[which(names(geno.dims) == "allele")]
  locus.dim <- geno.dims[which(names(geno.dims) == "locus")]
  
  
  ref.col <- which(dimnames(gene)[[allele.dim]] == ref.allele)
  new.allele.names <- dimnames(gene)[[allele.dim]][-ref.col]
  
  #=====================================================================
  #begin code for multi-allelic cross
  #=====================================================================
  
  #apply the modeling function to each marker column
  # !diagnostics suppress=m
  if (run.parallel) {
    
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    cape.dir.full <- find.package("cape")
    cape.dir <- str_replace(cape.dir.full,"cape_pkg/cape","cape_pkg")
    parallel::clusterExport(cl, "cape.dir", envir=environment())
    parallel::clusterEvalQ(cl, .libPaths(cape.dir))
    results <- foreach::foreach(m = 1:dim(gene)[[locus.dim]], .packages = 'cape', .export = c("get.stats.multiallele", "check.geno")) %dopar% {
      get.stats.multiallele(phenotype.vector, gene[,,m], covar.table = covar.table, 
      model.family, ref.col)
    }
    parallel::stopCluster(cl)
    
  } else {
    
    index <- 1:dim(gene)[[locus.dim]]
    results <- lapply(index, function(x) get.stats.multiallele(phenotype.vector, gene[,,x], 
      covar.table, model.family, ref.col))    
  }
  
  t.stat.array <- matrix(unlist(lapply(results, function(x) x[[1]]["t.stat",])), 
  ncol = length(new.allele.names), byrow = TRUE)
  colnames(t.stat.array) <- new.allele.names
  
  
  return(t.stat.array)		
}

