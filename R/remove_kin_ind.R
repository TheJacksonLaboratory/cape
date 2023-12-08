#' Removes individuals from the kinship object to match the cape.obj
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param kin_obj a kinship object
#' 
#' @return updated kinship object
#'
#' @export
remove_kin_ind <- function(data_obj, kin_obj){
  
  class_kin <- class(kin_obj)
  if(class_kin == "matrix"){
    ind_locale <- intersect(rownames(data_obj$pheno), rownames(kin_obj))
    new_kin_obj <- kin_obj[ind_locale, ind_locale]
  }else{
    new_kin_obj <- lapply(kin_obj, function(x) 
      x[intersect(rownames(data_obj$pheno), rownames(x)), intersect(rownames(data_obj$pheno),colnames(x))])
  }
  
  return(new_kin_obj)	
}