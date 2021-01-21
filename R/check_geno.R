#' Check to see if the any markers in the supplied
#' genotype matrix is one of the covariates.
#'
#' @param genotype The vector of genotype values for
#' a single marker
#' @param covar_table The table of covariates
#' 
#' @return covariates indices that match the genotype
#' @keywords internal
#'
check_geno <- function(genotype, covar_table){
  
  if(is.null(covar_table)){
    return(NULL)
  }
  
  covar_which <- NULL
  geno <- as.vector(genotype[,1])
  for(i in 1:ncol(covar_table)){
    
    if(identical(as.vector(covar_table[,i]), geno)){
      covar_which <- c(covar_which, i)
    }
  }
  
  return(covar_which)
}
