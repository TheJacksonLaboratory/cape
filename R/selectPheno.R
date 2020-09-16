#' This function selects the phenotypes in a Cape object
#' 
#' Updates the pheno object to include only `pheno.which` columns.
#' Optionally scale and/or normalize traits.
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param pheno.which vector of names from the parameters YAML file.
#' This vector should include both traits and covariates. The covariates
#' are assigned after trait selection. 
#' @param min.entries minimum number of data entries the phenotype needs 
#' to have for it to be included. If any trait has fewer than min.entries,
#' It will be removed with a warning.
#' @param scale.pheno if TRUE then phenotypes are mean-centered and standardized
#' @param rank.norm.pheno if TRUE then phenotypes are rank Z normalized
#'
#' @return updated \code{\link{Cape}} object
#' 
#' @export
selectPheno <- function(data.obj, pheno.which, min.entries = 5, scale.pheno = FALSE, rank.norm.pheno = FALSE){
  check.underscore(data.obj)
  # check.bad.markers(data.obj)
  
  
  pheno <- data.obj$pheno
  
  #find the phenotype column numbers if 
  #names have been put in instead of numbers	
  pheno.num <- get.col.num(pheno, pheno.which)
  
  if(length(pheno.num) < 2){
    stop("There must be at least two phenotypes selected.")
  }

  new.pheno <- pheno[,pheno.num]
  
  # #make sure the phenotypes are numeric
  # #and replace the phenotype matrix with 
  # #the selected phenotypes
  new.pheno <- apply(new.pheno, 2, as.numeric)
  
  #check to see if there are any phenotypes with
  #fewer than 5 entries
  data.entries <- as.vector(apply(new.pheno, 2, function(x) length(which(!is.na(x)))))
  bad.pheno <- which(data.entries <= min.entries)
  
  if(length(bad.pheno) > 0){
    final.pheno <- new.pheno[,-bad.pheno]
    cat("The following phenotypes had fewer than ", min.entries, " entries and were removed.\n")
    cat(paste("\t", colnames(new.pheno)[bad.pheno], "\n"))
  }else{
    final.pheno <- new.pheno
  }
  
  if(rank.norm.pheno){
    final.pheno <- apply(final.pheno, 2, rz.transform)
  }
  
  if(scale.pheno){
    final.pheno <- apply(final.pheno, 2, center.std) #mean center and standardize the phenotypes
  }
  
  rownames(final.pheno) <- rownames(pheno)		
  data.obj$pheno <- final.pheno
  
  return(data.obj)
  
}