#' Check to see if the any markers in the supplied
#' genotype matrix is one of the covariates.
#'
#' @param genotype The vector of genotype values for
#' a single marker
#' @param covar.table The table of covariates
#' 
#' @return covariates indices that match the genotype
#'
check.geno <- function(genotype, covar.table){
  
  if(is.null(covar.table)){
    return(NULL)
  }
  
  covar.which <- NULL
  geno <- as.vector(genotype[,1])
  for(i in 1:ncol(covar.table)){
    
    if(identical(as.vector(covar.table[,i]), geno)){
      covar.which <- c(covar.which, i)
    }
  }
  
  return(covar.which)
}
