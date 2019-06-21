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
one.singlescanDO <- function(phenotype.vector, genotype.mat, model.family, ref.allele = "A", covar.vector = NULL, run.parallel = TRUE, n.cores = 4){
  
  if(!run.parallel){n.cores = 1}
  
  gene <- genotype.mat
  
  #Get the dimension names to minimize confusion	
  mouse.dim <- which(names(dimnames(gene)) == "mouse")
  locus.dim <- which(names(dimnames(gene)) == "locus")
  allele.dim <- which(names(dimnames(gene)) == "allele")
  
  
  #if there are covariates specified, pull these out.
  #covariates must be coded as markers and contained in the
  #genotype matrix
  if(!is.null(covar.vector)){
    covar <- names(covar.vector[which(covar.vector == 1)])
    if(length(covar) > 0){
      covar.loc <- get.col.num(gene, covar, locus.dim)
      covar.table <- array(NA, dim = c(dim(gene)[[mouse.dim]], dim(gene)[[allele.dim]], length(covar)))
      covar.table[,,1:length(covar)] <- gene[,,covar.loc]
    }else{
      covar.table = NULL
    }
  }else{
    covar.table <- NULL
  }
  
  ref.col <- which(dimnames(gene)[[allele.dim]] == ref.allele)
  new.allele.names <- dimnames(gene)[[allele.dim]][-ref.col]
  
  #=====================================================================
  #begin code for multi-allelic cross
  #=====================================================================
  
  #apply the modeling function to each marker column
  if (run.parallel) {
    
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    results <- foreach::foreach(m = 1:dim(gene)[[locus.dim]], .export = c("get.stats.multiallele", "check.geno")) %dopar% {
      get.stats.multiallele(phenotype.vector, gene[,,m], covar.table, model.family, ref.col)
    }
    parallel::stopCluster(cl)
    
  } else {
    
    results <- c()
    index <- 1:dim(gene)[[locus.dim]]
    for (m in index) {
      results <- cbind(results, get.stats.multiallele(phenotype.vector, gene[,,m], covar.table, model.family, ref.col))
    }
    
  }
  
  t.stat.array <- matrix(unlist(lapply(results, function(x) x[[1]]["t.stat",])), ncol = length(new.allele.names), byrow = TRUE)
  
  colnames(t.stat.array) <- new.allele.names
  
  
  return(t.stat.array)		
}

