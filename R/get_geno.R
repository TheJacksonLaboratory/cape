#' Gets the geno object
#'
#' This is an internal function returns the 
#' genotype matrix for scanning as defined by 
#' the markers and individuals specified in 
# data_obj$geno_names.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object.
#'
#' @return Returns the genotype array matching 
#' the markers and individuals specified in 
#' data_obj$geno_names
#'
#' @importFrom abind abind
#' 
#' @examples 
#' \dontrun{
#' geno <- get_geno(data_obj, geno_obj)
#' }
#'
#' @export

get_geno <- function(data_obj, geno_obj){
    
  geno_names <- data_obj$geno_names
  
  if(is.null(geno_obj)){
    geno <- data_obj$geno
  }else{
    if(class(geno_obj) == "array"){
      geno <- geno_obj
    }else{
      geno <- geno_obj$geno	
    }
  }
  
  
  if(is.null(geno)){
    stop("I can't find the genotype data. Please make sure it is in either data_obj or geno_obj.")
  }
  
  geno_dims <- get_geno_dim()
  mouse_dim <- geno_dims[which(names(geno_dims) == "mouse")]
  allele_dim <- geno_dims[which(names(geno_dims) == "allele")]
  locus_dim <- geno_dims[which(names(geno_dims) == "locus")]
  
  #subset the genotype object to match the 
  #individuals and markers we want to scan
  ind_locale <- match(geno_names[[mouse_dim]], dimnames(geno)[[mouse_dim]])
  allele_locale <- match(geno_names[[allele_dim]], dimnames(geno)[[allele_dim]])
  locus_locale <- match(geno_names[[locus_dim]], dimnames(geno)[[locus_dim]])
  
  #check for NAs, meaning the locus from the data object cannot be
  #found in the genotyope object
  locus_locale <- locus_locale[which(!is.na(locus_locale))]
  
  gene <- geno[ind_locale, allele_locale, locus_locale]
  
  #if there is a covariate table in the data object, this is added
  #to the genotype object
  
  if(!is.null(data_obj$covar_table)){
    covar_vals <- data_obj$covar_table
    covar_names <- colnames(covar_vals)
    covar_table <- array(NA, dim = c(length(geno_names[[mouse_dim]]), length(geno_names[[allele_dim]]), dim(covar_vals)[2]))
    for(i in 1:dim(covar_vals)[2]){
      covar_table[,1:dim(covar_table)[2],i] <- covar_vals[,i]
    }
    dimnames(covar_table)[[3]]  <- covar_names
    gene <- abind(gene, covar_table, along = 3)
  }
  
  names(dimnames(gene)) <- names(dimnames(geno))
  
  return(gene)	
}