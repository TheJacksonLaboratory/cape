#' Return the genotype matrix with covariates 
#' added. 
#' 
#' For pairwise testing, cape appends the covariates
#' to the genotype matrix. This function performs
#' that task.
#' 
#' @param data_obj The cape object. See \code{\link{Cape}}.
#' @param geno_obj A genotype object.
#' @param g_covar A logical value indicating whether to append
#' genotype-derived covariates to the final matrix. Defaults to TRUE.
#' @param p_covar A logical value indicating whether to append
#' phenotype-derived covariates to the final matrix. Defaults to TRUE.
#' @param for_pairscan If TRUE, the function returns the genotype matrix 
#' already designated for the pairscan. Otherwise it returns the full
#' genotype matrix as specified by data_obj$geno_names.
#'
#' @return Returns a genotype matrix with covariates appended.
#' @keywords internal


#This function puts covariates into the genotype
#matrix for easier testing of all pairs

get_geno_with_covar <- function(data_obj, geno_obj = NULL, g_covar = TRUE, p_covar = TRUE, for_pairscan = TRUE){
  
  
  covar_info <- get_covar(data_obj)
  
  covar_locale <- NULL
  if(g_covar){
    covar_locale <- c(covar_locale, which(covar_info$covar_type == "g"))
  }
  if(p_covar){
    covar_locale <- c(covar_locale, which(covar_info$covar_type == "p"))
  }
  
  
  if(for_pairscan){
    geno <- data_obj$geno_for_pairscan			
  }else{
    geno <- get_geno(data_obj, geno_obj)
  }
  
  is_char <- as.logical(is.na(suppressWarnings(as.numeric(colnames(geno)[1]))))
  
  if(is_char){
    colnames(covar_info$covar_table) <- covar_info$covar_names
  }
  
  geno <- cbind(geno, covar_info$covar_table[,covar_locale,drop=FALSE])
  
  #if there are marker covariates make sure these
  #are placed in the right order in the genotype
  #matrix
  if(g_covar && length(which(covar_info$covar_type == "g"))  > 0){
    new_geno_chr <- get_marker_chr(data_obj, colnames(geno))
    new_geno_pos <- get_marker_location(data_obj, colnames(geno))
    marker_pos_table <- cbind(new_geno_chr, new_geno_pos)
    marker_pos_table <- sort_by_then_by(marker_pos_table, col_type = c("n", "n"), return_order = TRUE)
    for(i in 1:2){
      geno <- geno[,marker_pos_table[,i]]
    }
    
    #now we need to make sure any phenotypic covariates
    #are at the end
    new_geno_chr <- get_marker_chr(data_obj, colnames(geno))
    pheno_covar_locale <- which(new_geno_chr == 0)
    if(length(pheno_covar_locale) > 0){
      geno <- cbind(geno[,-pheno_covar_locale,drop=FALSE],geno[,pheno_covar_locale,drop=FALSE])
    }
  }	
  
  return(geno)
  
}


