#' Create a covariate from a trait
#' 
#' This function takes a variable from the phenotype matrix
#' for example, diet treament or sex and converts it to
#' a covariate.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param pheno_which vector of trait names to be used as covariates
#'
#' @return Returns the data object with the specified traits removed
#' from the phenotype matrix and transferred where they will be used
#' as covariates. Information about assigned covariates can be retrieved
#' with \code{\link{get_covar}}.
#'
#' @export

pheno2covar <- function(data_obj, pheno_which){
  
  pheno <- data_obj$pheno
  if(length(rownames(pheno)) == 0){stop("The phenotype matrix must have rownames.")}
  
  marker_locale <- get_col_num(pheno, pheno_which)
  
  if(length(marker_locale) == 0){
    return(data_obj)
  }
  
  if(length(unique(marker_locale)) != length(pheno_which)){
    stop("Phenotype labels cannot be substrings of other phenotypes.")
  }
  
  
  #make a separate covariate table, then modify the dimnames
  #in the genotype object to include the covariates
  #do not modify the genotype object
  
  #get information for any covariates that have already been specified
  covar_info <- get_covar(data_obj)
  
  covar_table <- pheno[,marker_locale,drop=FALSE]
  rownames(covar_table) <- rownames(pheno)
  
  #check for invariant covariates
  covar_var <- apply(covar_table, 2, function(x) var(x, na.rm = TRUE))
  if(any(covar_var == 0)){
    cat("Removing invariant covariates...\n")
    covar_table <- covar_table[,-which(covar_var == 0),drop=FALSE]
    data_obj$covariates <- data_obj$covariates[-which(covar_var == 0)]
  }

  #if after removing invariant covariates, there
  #are none left, just return the data object.
  if(length(data_obj$covariates) == 0){
    return(data_obj)
  }

  #the covariate's number starts after genetic markers and any existing phenotypic covariates
  start_covar <- max(as.numeric(data_obj$marker_num))+1+length(covar_info$covar_names)
  colnames(covar_table) <- start_covar:(start_covar+dim(covar_table)[2]-1)
  
  #scale the covariate(s) to be between 0 and 1
  scaled_table <- apply(covar_table, 2, function(x) x + abs(min(x, na.rm = TRUE)))
  scaled_table <- apply(scaled_table, 2, function(x) x/max(x, na.rm = TRUE))

  #add the covariates to any existing covariates
  data_obj$p_covar_table <- cbind(data_obj$p_covar_table, scaled_table)
  
  #take the phenotypes made into markers out of the phenotype matrix
  new_pheno <- pheno[,-marker_locale,drop=FALSE]
  data_obj$pheno <- new_pheno
  
  if (is.null(data_obj$p_covar)) {
    data_obj$p_covar <- data_obj$covariates
  } else {
    data_obj$p_covar <- c(data_obj$p_covar, pheno_which)
  }
  
  return(data_obj)
  
}