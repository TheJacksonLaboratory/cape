#' Add eigentraits, singular values, and right singular vectors to the Cape object
#' 
#' This script takes in the data object and returns 
#' the it with eigentraits added to the data object
#' along with the singular values and the right
#' singular values.
#' if scale.pheno is TRUE, the phenotypes are mean
#' centered and standardized before the svd is run
#' The user also has the option to rank normalize 
#' the phenotypes. This argument defaults to false
#' because many users will do their own normalization.
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param scale.pheno boolean, default: TRUE
#' @param normalize.pheno boolean, default: TRUE
#'
#' @return updated \code{\link{Cape}} object
#'
#' @export
get.eigentraits <- function(data.obj, scale.pheno = TRUE, normalize.pheno = TRUE){
  
  
  if(ncol(data.obj$pheno) <= 1){
    stop("CAPE requires more than one phenotype.")
  }
  
  #first make sure there are no individuals
  #with missing phenotypes. This also makes 
  #sure the phenotypes are numeric
  ind.missing.pheno <- which(is.na(data.obj$pheno), arr.ind = TRUE)
  if(nrow(ind.missing.pheno) > 0){
    cat("Removing ", length(unique(ind.missing.pheno[,1])), " individuals with missing phenotypes.")
  }
  data.obj <- remove.ind(data.obj, ind.to.remove = unique(ind.missing.pheno[,1]))
  
  #perform a variance check on the new covariates
  
  if(!is.null(data.obj$p_covar_table)){
    var.check <- apply(data.obj$p_covar_table, 2, function(x) var(x, na.rm = TRUE))
    if(any(var.check == 0)){
      zero.locale <- which(var.check == 0)
      cat("Some covariates now have zero variance. Removing:")
      cat(data.obj$p_covar[zero.locale], sep = "\n")
      data.obj$p_covar_table <- as.array(data.obj$p_covar_table[,which(var.check > 0)])
      data.obj$p_covar <- as.array(data.obj$p_covar[which(var.check > 0)])
    }
    
    # TODO this code is duplicated in singlescan
    
    #also remove the NAs and check the matrix for rank
    not.na.locale <- which(!is.na(apply(data.obj$p_covar_table,1,sum)))
    no.na.cov <- as.array(data.obj$p_covar_table[not.na.locale,drop=FALSE])
    design.cov <- cbind(rep(1, dim(no.na.cov)[1]), no.na.cov)  #####################
    rank.cov <- Matrix::rankMatrix(design.cov)
    if(rank.cov[[1]] < dim(design.cov)[2]){
      stop("The covariate matrix does not appear to be linearly independent.\nIf you are using dummy variables for groups, leave one of the groups out.")
    }
  }
  
  
  pheno <- data.obj$pheno
  
  if(dim(pheno)[1] < dim(pheno)[2]){ #if there are more individuals than phenotypes
    stop("The system is underdetermined. There are more phenotypes than individuals.")
  }
  
  #This function mean centers and standardizes a vector
  center.std <- function(v){
    mean.v <- mean(v)
    centered <- v - mean.v
    sd.v <- sd(v)
    final.v <- centered/sd.v
    return(final.v)
  }
  
  if(normalize.pheno){
    pheno <- apply(pheno, 2, rz.transform)
  }
  
  if(scale.pheno){
    pheno <- apply(pheno, 2, center.std) #mean center and standardize the phenotypes
  }
  
  
  data.obj$pheno <- pheno #replace the raw phenotypes with scaled, normalized phenotypes (if we have done those things)
  
  svd.pheno <- svd(pheno)
  
  
  #add the eigentraits and singular values to the data object
  
  data.obj$ET <- svd.pheno$u
  rownames(data.obj$ET) <- rownames(data.obj$pheno)
  data.obj$right_singular_vectors <- svd.pheno$v
  data.obj$singular_values <- svd.pheno$d
  data.obj$traits_scaled <- as.logical(scale.pheno)
  data.obj$traits_normalized <- as.logical(normalize.pheno)
  
  
  
  colnames(data.obj$ET) <- paste("ET", 1:length(svd.pheno$u[1,]), sep = "")
  
  return(data.obj)
  
}