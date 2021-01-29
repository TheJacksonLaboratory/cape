#' This function selects the phenotypes in a Cape object
#' 
#' Updates the pheno object to include only `pheno_which` columns.
#' Optionally scale and/or normalize traits.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param pheno_which vector of names from the parameters YAML file.
#' This vector should include both traits and covariates. The covariates
#' are assigned after trait selection. 
#' @param min_entries minimum number of data entries the phenotype needs 
#' to have for it to be included. If any trait has fewer than min_entries,
#' It will be removed with a warning.
#' @param scale_pheno if TRUE then phenotypes are mean-centered and standardized
#' @param rank_norm_pheno if TRUE then phenotypes are rank Z normalized
#'
#' @return updated \code{\link{Cape}} object
#' 
#' @examples 
#' \dontrun{
#' data_obj <- select_pheno(data_obj, pheno_which = c("BW_24", "INS_24", "log_GLU_24"))
#' }
#' 
#' @export
select_pheno <- function(data_obj, pheno_which, min_entries = 5, scale_pheno = FALSE, rank_norm_pheno = FALSE){
  check_underscore(data_obj)
  # check_bad_markers(data_obj)
  
  
  pheno <- data_obj$pheno
  
  #find the phenotype column numbers if 
  #names have been put in instead of numbers	
  pheno_num <- get_col_num(pheno, pheno_which)
  
  if(length(pheno_num) < 2){
    stop("There must be at least two phenotypes selected.")
  }

  new_pheno <- pheno[,pheno_num]
  
  # #make sure the phenotypes are numeric
  # #and replace the phenotype matrix with 
  # #the selected phenotypes
  new_pheno <- apply(new_pheno, 2, as.numeric)
  
  #check to see if there are any phenotypes with
  #fewer than 5 entries
  data_entries <- as.vector(apply(new_pheno, 2, function(x) length(which(!is.na(x)))))
  bad_pheno <- which(data_entries <= min_entries)
  
  if(length(bad_pheno) > 0){
    final_pheno <- new_pheno[,-bad_pheno]
    message("The following phenotypes had fewer than ", min_entries, " entries and were removed:\n", paste(colnames(new_pheno)[bad_pheno], collapse = ", "))
  }else{
    final_pheno <- new_pheno
  }
  
  if(rank_norm_pheno){
    final_pheno <- apply(final_pheno, 2, rz_transform)
  }
  
  if(scale_pheno){
    final_pheno <- apply(final_pheno, 2, center_std) #mean center and standardize the phenotypes
  }
  
  rownames(final_pheno) <- rownames(pheno)		
  data_obj$pheno <- final_pheno
  
  return(data_obj)
  
}