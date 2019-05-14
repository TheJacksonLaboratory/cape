#' Check to see if the genotype matrix is one of the covariates
#'
#' @param genotype
#' @param covar.table
#' 
#' @return covariates indices that match the genotype
#'
#' @export
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
