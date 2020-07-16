#' Corrects genotypes, phenotypes, and covariates
#' for kinship.
#' 
#' This function uses linear mixed models to adjust the
#' genotype matrix, phenotype matrix, and covariate matrix
#' for kinship based on the kinship matrix calculated 
#' by \link{\code{Kinship}}.
#'
#' @param kin.obj The kinship object calculated by \link{\code{Kinship}}
#' @param geno a genotype object.
#' @param chr1 The first of two chromomsomes to leave out of the calculation, if any.
#' @param chr2 The second of two chromomsomes to leave out of the calculation, if any.
#' @param phenoV The phenotype vector
#' @param covarV The covariate vector or matrix
#'
#' @details If using leave-one-chromosome-out (LOCO), chr1 and chr2 should have the same value.
#' If chr1 and chr2 are different, both specified, leave-two-chromosomes-out (LTCO) will be
#' used. After testing LTCO, we do not recommend using this method. We instead recommend 
#' always using the overall kinship correction. In this case, both chr1 and chr2 should
#' be NULL, their default values.

#' @return This function returns a list with the corrected phenotype,
#' genotypes, and covariates. These are used in \link{\code{singlescan}}
#' and \link{\code{pairscan.kin}}.
#'

kinship.on.the.fly <- function(kin.obj, geno, chr1 = NULL, chr2 = NULL, phenoV = NULL, 
covarV = NULL, verbose = FALSE){
  
  get.g=function(pair = NULL, phenotype, covarV){

    if(is.null(pair) || unique(pair) == "overall"){pair = NULL}
    
    if(is.null(pair)){
      pair.name <- "overall"
    }else{
      pair.name <- paste(pair, collapse = ",")
    }
    
    if(verbose){cat("Chromosomes:", pair.name, "\n")}
    
    if(class(kin.obj) == "matrix"){
      full.kin <- kin.obj
    }else{
      kin.mat.locale <- which(names(kin.obj) == pair.name)
      full.kin <- kin.obj[[kin.mat.locale]]
    }
    
    #remove individuals with NAs
    is.na.pheno <- which(is.na(phenotype))
    if(length(covarV) > 0){
      is.na.covar <- unique(which(is.na(covarV), arr.ind = TRUE)[,1])
      }else{
        is.na.covar <- NULL
      }
    all.na <- unique(c(is.na.pheno, is.na.covar))
    not.na <- setdiff(1:length(phenotype), all.na)
    no.na.ind <- rownames(phenotype)[not.na]
    common.ind <- intersect(no.na.ind, colnames(full.kin))

    kin.locale <- match(common.ind, colnames(full.kin))  
    K <- full.kin[kin.locale,kin.locale]
    pheno.locale <- match(common.ind, rownames(phenotype))

    #for the corrections below, look into including epistatic kinship 
    #matrices. This may help us gain power to see epistatic interactions

    #if we are correcting the covariate only don't put it in the model
    if(verbose){cat("\tFitting model...\n")}
    if(is.null(covarV) || is.null(pair)){
      model = regress(as.vector(phenotype[pheno.locale])~1,~K, pos = c(TRUE, TRUE), 
      tol = 1e-2)
    }else{
      model = regress(as.vector(phenotype)[pheno.locale]~covarV[pheno.locale,], ~K, 
      pos = c(TRUE, TRUE), tol = 1e-2)
    }
    
    #This err.cov is the same as err.cov in Dan's code using estVC
    #err.cov = summary(model)$sigma[1]*K+summary(model)$sigma[2]*diag(nrow(K))
    if(verbose){cat("\tCalculating err.cov...\n")}
    err.cov = model$sigma[1]*K+model$sigma[2]*diag(nrow(K))
    
  if(verbose){cat("\tCalculating eW...\n")}
    eW = eigen(err.cov, symmetric = TRUE)
    if(min(eW$values) < 0 && abs(min(eW$values)) > sqrt(.Machine$double.eps)){
    }else{
      eW$values[eW$values <= 0] = Inf
    } 
    err.cov = eW$vector %*% diag(eW$values^-0.5) %*% t(eW$vector)
    
    new.pheno <- err.cov %*% phenotype[pheno.locale,]
    
    if(length(dim(geno)) == 3){
      l.geno <- lapply(1:dim(geno)[2], function(x) err.cov %*% geno[pheno.locale,x,]); #hist(new.geno)
      new.geno <- array(NA, dim = dim(geno[pheno.locale,,]))
      for(i in 1:length(l.geno)){
        new.geno[,i,] <- l.geno[[i]]
      }
      dimnames(new.geno) <- dimnames(geno[pheno.locale,,])
    }else{
      new.geno <- err.cov %*% geno[pheno.locale,]
      dimnames(new.geno) <- dimnames(geno[pheno.locale,])				
    }
    
    if(!is.null(covarV)){
      new.covar <- err.cov %*% covarV[pheno.locale,]
    }else{
      new.covar <- NULL	
    }
    
    results = list(err.cov, new.pheno, new.geno, new.covar)
    names(results) <- c("err.cov", "corrected.pheno", "corrected.geno", "corrected.covar")
    return(results)
  }
  

  result <- get.g(pair = c(chr1, chr2), phenotype = phenoV, covarV = covarV)
  
  
  return(result)
  
}
