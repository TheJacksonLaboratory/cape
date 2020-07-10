#' Performs a 1D scan for generating the 2D null distribution
#' 
#' @param phenotype.vector A vector of phenotype values, one entry for each individual.
#' @param genotype.mat A matrix of genotype values with individuals in rows and markers 
#' in columns. Matrix entries contain the probability of the reference allele at each 
#' position for each individual.
#' @param model.family Indicates the model family of the phenotypes. This can be 
#'   either "gaussian" or "binomial".
#' @param ref.allele the reference allele from the Cape data.obj
#' @param covar.vector A vector of 1's and 0's indicating which markers should be used 
#' as covariats (1) and which should not (0), optional
#' @param run.parallel default = TRUE
#' @param n.cores default = 4
#' 
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

