#' Removes individuals from the kinship object to match the cape.obj
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param kin.obj a kinship object
#'
#' @return updated kinship object
#'
#' @export
remove.kin.ind <- function(data.obj, kin.obj){
  
  if(class(kin.obj) == "matrix"){
    ind.locale <- intersect(rownames(data.obj$pheno), rownames(kin.obj))
    new.kin.obj <- kin.obj[ind.locale, ind.locale]
  }else{
    new.kin.obj <- lapply(kin.obj, function(x) 
      x[intersect(rownames(data.obj$pheno), rownames(x)), intersect(rownames(data.obj$pheno),colnames(x))])
  }
  
  return(new.kin.obj)	
}