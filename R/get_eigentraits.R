#' Calculate eigentraits
#' 
#' This function uses singular value decomposition
#' (SVD) to calculate eigentraits from the phenotype 
#' matrix in the cape data object. It
#' adds the eigentrait matrix to the data object
#' along with the singular values and the right
#' singular vectors.
#' 
#' If scale_pheno is TRUE, the phenotypes are 
#' mean-centered and standardized before running
#' the svd.
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param scale_pheno A logical value indicating whether to 
#' mean-center and standardize the traits.
#' @param normalize_pheno A logical value indicating whether to 
#' rankZ normalize the phenotypes.
#'
#' @details Because we use SVD in this step, there can be
#' no missing values in the phenotype matrix. Any individuals
#' with missing values are removed with a message. 

#' @return Returns the data object with the eigentraits,
#' singular values, and right singular vectors added.
#'
#' @importFrom Matrix rankMatrix
#'
#' @export

get_eigentraits <- function(data_obj, scale_pheno = TRUE, normalize_pheno = TRUE){
  
  
  if(ncol(data_obj$pheno) <= 1){
    stop("CAPE requires more than one phenotype.")
  }
  
  #first make sure there are no individuals
  #with missing phenotypes. This also makes 
  #sure the phenotypes are numeric
  ind_missing_pheno <- which(is.na(data_obj$pheno), arr.ind = TRUE)
  if(nrow(ind_missing_pheno) > 0){
    message("Removing ", length(unique(ind_missing_pheno[,1])), " individuals with missing phenotypes.\n")
  }
  data_obj <- remove_ind(data_obj, ind_to_remove = unique(ind_missing_pheno[,1]))
  
  #perform a variance check on the new covariates
  
  if(!is.null(data_obj$p_covar_table)){
    var_check <- apply(data_obj$p_covar_table, 2, function(x) var(x, na.rm = TRUE))
    if(any(var_check == 0)){
      zero_locale <- which(var_check == 0)
      message("Some covariates now have zero variance.\nRemoving: ", paste(data_obj$p_covar[zero_locale], collapse = ", "))
      data_obj$p_covar_table <- as.array(data_obj$p_covar_table[,which(var_check > 0)])
      data_obj$p_covar <- as.array(data_obj$p_covar[which(var_check > 0)])
    }
    
    # TODO this code is duplicated in singlescan
    
    #also remove the NAs and check the matrix for rank
    not_na_locale <- which(!is.na(apply(data_obj$p_covar_table,1,sum)))
    no_na_cov <- as.array(data_obj$p_covar_table[not_na_locale,drop=FALSE])
    design_cov <- cbind(rep(1, dim(no_na_cov)[1]), no_na_cov)  #####################
    rank_cov <- rankMatrix(design_cov)
    if(rank_cov[[1]] < dim(design_cov)[2]){
      stop("The covariate matrix does not appear to be linearly independent.\nIf you are using dummy variables for groups, leave one of the groups out.")
    }
  }
  
  
  pheno <- data_obj$pheno
  
  if(dim(pheno)[1] < dim(pheno)[2]){ #if there are more individuals than phenotypes
    stop("The system is underdetermined. There are more phenotypes than individuals.")
  }
  
  #This function mean centers and standardizes a vector
  center_std <- function(v){
    mean_v <- mean(v)
    centered <- v - mean_v
    sd_v <- sd(v)
    final_v <- centered/sd_v
    return(final_v)
  }
  
  #preserve the initial phentype matrix.
  data_obj$raw_pheno <- pheno 

  if(normalize_pheno){
    pheno <- apply(pheno, 2, rz_transform)
  }
  
  if(scale_pheno){
    pheno <- apply(pheno, 2, center_std) #mean center and standardize the phenotypes
  }
  
  data_obj$pheno <- pheno #replace the raw phenotypes with scaled, normalized phenotypes (if we have done those things)
  
  svd_pheno <- svd(pheno)
  
  #add the eigentraits and singular values to the data object
  
  data_obj$ET <- svd_pheno$u
  rownames(data_obj$ET) <- rownames(data_obj$pheno)
  data_obj$right_singular_vectors <- svd_pheno$v
  data_obj$singular_values <- svd_pheno$d
  data_obj$traits_scaled <- as.logical(scale_pheno)
  data_obj$traits_normalized <- as.logical(normalize_pheno)
  
  
  
  colnames(data_obj$ET) <- paste("ET", 1:length(svd_pheno$u[1,]), sep = "")
  
  return(data_obj)
  
}