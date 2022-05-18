#' Corrects genotypes, phenotypes, and covariates
#' for kinship.
#' 
#' This function uses linear mixed models to adjust the
#' genotype matrix, phenotype matrix, and covariate matrix
#' for kinship based on the kinship matrix calculated 
#' by \code{\link{kinship}}.
#'
#' @param kin_obj The kinship object calculated by \code{\link{kinship}}
#' @param geno a genotype object.
#' @param chr1 The first of two chromosomes to leave out of the calculation, if any.
#' @param chr2 The second of two chromosomes to leave out of the calculation, if any.
#' @param phenoV The phenotype vector
#' @param covarV The covariate vector or matrix
#' @param verbose A logical value indicating whether to print progress to the screen
#'
#' @details If using leave-one-chromosome-out (LOCO), chr1 and chr2 should have the same value.
#' If chr1 and chr2 are different, both specified, leave-two-chromosomes-out (LTCO) will be
#' used. After testing LTCO, we do not recommend using this method. We instead recommend 
#' always using the overall kinship correction. In this case, both chr1 and chr2 should
#' be NULL, their default values.

#' @return This function returns a list with the corrected phenotype,
#' genotypes, and covariates. These are used in \code{\link{singlescan}}
#' and \code{\link{pairscan_kin}}.
#'
#' @importFrom regress regress
#' @keywords internal

kin_adjust <- function(kin_obj, geno, chr1 = NULL, chr2 = NULL, phenoV = NULL, 
covarV = NULL, verbose = FALSE){
  
  get_g=function(pair = NULL, phenotype, covarV){

    if(is.null(pair) || unique(pair) == "overall"){pair = NULL}
    
    if(is.null(pair)){
      pair_name <- "overall"
    }else{
      pair_name <- paste(pair, collapse = ",")
    }
    
    if(verbose){cat("Chromosomes:", pair_name, "\n")}
    
    class_kin <- class(kin_obj)[1]
    if(class_kin == "matrix"){
      full_kin <- kin_obj
      if(verbose){cat("\tUsing overall matrix\n")}
    }else{
      kin_mat_locale <- which(names(kin_obj) == pair_name)
      full_kin <- kin_obj[[kin_mat_locale]]
      if(verbose){cat("\tUsing", pair_name, "kinship matrix\n")}
    }
    
    #remove individuals with NAs
    is_na_pheno <- which(is.na(phenotype))
    if(length(covarV) > 0){
      is_na_covar <- unique(which(is.na(covarV), arr.ind = TRUE)[,1])
      }else{
        is_na_covar <- NULL
      }
    all_na <- unique(c(is_na_pheno, is_na_covar))
    not_na <- setdiff(1:length(phenotype), all_na)
    no_na_ind <- rownames(phenotype)[not_na]
    common_ind <- intersect(no_na_ind, colnames(full_kin))

    kin_locale <- match(common_ind, colnames(full_kin))  
    K <- full_kin[kin_locale,kin_locale]
    pheno_locale <- match(common_ind, rownames(phenotype))

    #for the corrections below, look into including epistatic kinship 
    #matrices. This may help us gain power to see epistatic interactions

    #if we are correcting the covariate only don't put it in the model
    if(verbose){cat("\tFitting model...\n")}
    if(is.null(covarV) || is.null(pair)){
      model = regress(as.vector(phenotype[pheno_locale])~1,~K, pos = c(TRUE, TRUE), 
      tol = 1e-2)
    }else{
      model = regress(as.vector(phenotype)[pheno_locale]~covarV[pheno_locale,], ~K, 
      pos = c(TRUE, TRUE), tol = 1e-2)
    }
    
    #This err_cov is the same as err_cov in Dan's code using estVC
    #err_cov = summary(model)$sigma[1]*K+summary(model)$sigma[2]*diag(nrow(K))
    #if(verbose){cat("\tCalculating err_cov...\n")}
    err_cov = model$sigma[1]*K+model$sigma[2]*diag(nrow(K))
    
  if(verbose){cat("\tCalculating eW...\n")}
    eW = eigen(err_cov, symmetric = TRUE)
    if(min(eW$values) < 0 && abs(min(eW$values)) > sqrt(.Machine$double.eps)){
    }else{
      eW$values[eW$values <= 0] = Inf
    } 
    err_cov = eW$vector %*% diag(eW$values^-0.5) %*% t(eW$vector)
    
    new_pheno <- err_cov %*% phenotype[pheno_locale,]
    
    if(length(dim(geno)) == 3){
      l_geno <- lapply(1:dim(geno)[2], function(x) err_cov %*% geno[pheno_locale,x,]); #hist(new_geno)
      new_geno <- array(NA, dim = dim(geno[pheno_locale,,]))
      for(i in 1:length(l_geno)){
        new_geno[,i,] <- l_geno[[i]]
      }
      dimnames(new_geno) <- dimnames(geno[pheno_locale,,])
    }else{
      new_geno <- err_cov %*% geno[pheno_locale,]
      dimnames(new_geno) <- dimnames(geno[pheno_locale,])				
    }
    
    if(!is.null(covarV)){
      new_covar <- err_cov %*% covarV[pheno_locale,]
    }else{
      new_covar <- NULL	
    }
    
    results = list(err_cov, new_pheno, new_geno, new_covar)
    names(results) <- c("err_cov", "corrected_pheno", "corrected_geno", "corrected_covar")
    return(results)
  }
  
  if(verbose){cat("Chromosomes", chr1, chr2, "\n")}

  chr_pair <- c(chr1, chr2)
  result <- get_g(pair = chr_pair, phenotype = phenoV, covarV = covarV)
  
  
  return(result)
  
}
