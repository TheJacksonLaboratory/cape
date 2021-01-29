#' Create a covariate from a trait
#' 
#' This function takes a variable from the phenotype matrix
#' for example, diet treatment or sex and converts it to
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
#' @examples 
#' \dontrun{
#' #convert weight to a covariate
#' data_obj <- pheno2covar(data_obj, "weight")
#' }
#' @export

pheno2covar <- function(data_obj, pheno_which){
  
  pheno <- data_obj$pheno
  p_covar <- pheno_which
  if(length(rownames(pheno)) == 0){stop("The phenotype matrix must have rownames.")}
  
  marker_locale <- get_col_num(pheno, pheno_which)
  
  if(length(marker_locale) == 0){
    return(data_obj)
  }
  
  if(length(unique(marker_locale)) != length(pheno_which)){
    stop("Phenotype labels cannot be substrings of other phenotypes.")
  }
  
  #make a separate covariate table, then modify the dimnames
  #in the genotype object to include the covariates.
  #Do not modify the genotype object
  
  #get information for any covariates that have already been specified
  covar_info <- get_covar(data_obj)
  
  covar_table <- pheno[,marker_locale,drop=FALSE]
  rownames(covar_table) <- rownames(pheno)
  
  #check for invariant covariates
  covar_var <- apply(covar_table, 2, function(x) var(x, na.rm = TRUE))
                     
  if(all(is.na(covar_var))){
  	warning("The covariates are invariant.")
  	return(data_obj)
  }
  if(all(covar_var == 0)){
    warning("The covariates are invariant.")    
    return(data_obj)
  }

  if(any(covar_var == 0)){
    message("Removing invariant covariates...")
    invar.locale <- which(covar_var == 0)
    covar_table <- covar_table[,-invar.locale,drop=FALSE]
    p_covar <- p_covar[-invar.locale]
  }


  #the covariate's number starts after genetic markers and any existing phenotypic covariates
  start_covar <- max(as.numeric(data_obj$marker_num))+1+length(covar_info$covar_names)
  colnames(covar_table) <- start_covar:(start_covar+dim(covar_table)[2]-1)
  
  #scale the covariate(s) to be between 0 and 1
  scaled_table <- apply(covar_table, 2, function(x) x + abs(min(x, na.rm = TRUE)))
  scaled_table <- apply(scaled_table, 2, function(x) x/max(x, na.rm = TRUE))

  #add the covariates to the data object
  data_obj$p_covar_table <- scaled_table
  data_obj$p_covar <- p_covar

  #take the phenotypes made into covariates out of the phenotype matrix
  new_pheno <- pheno[,-marker_locale,drop=FALSE]
  data_obj$pheno <- new_pheno
    
  return(data_obj)
  
}