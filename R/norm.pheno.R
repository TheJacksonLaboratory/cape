#This function mean centers and rank z normalizes phenotypes
norm.pheno <- function(data.obj, mean.center = TRUE){
  
  pheno <- data.obj$pheno
  raw.pheno <- data.obj$pheno #retain the raw phenotype values
  
  pheno <- apply(pheno, 2, rz.transform)
  
  if(mean.center){
    pheno <- apply(pheno, 2, center.std) #mean center and standardize the phenotypes
  }
  
  data.obj$pheno <- pheno #replace the raw phenotypes with scaled, normalized phenotypes (if we have done those things)
  data.obj$raw_pheno <- raw.pheno
  
  return(data.obj)
}