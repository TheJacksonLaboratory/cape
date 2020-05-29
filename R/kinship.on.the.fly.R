#' Calculates the leave-two-out kinship matrix
#' 
#' This function calculates the leave-two-out kinship matrix
#' used in capeRel
#' plot.adj.mat is used to plot the covariance matrices 
#' vs. the positive definite matrices
#' to correct just for the covariate, leave chr1 and chr2 NULL
#'
#' @param kin.obj a kinship object
#' @param geno a genotype object. If this is not supplied then it is generated here.
#' @param chr1
#' @param chr2
#' @param phenoV
#' @param covarV
#'
#' @return \code{list("err.cov", "corrected.pheno", "corrected.geno", "corrected.covar")}
#'
#' @export
kinship.on.the.fly <- function(kin.obj, geno, chr1 = NULL, chr2 = NULL, phenoV = NULL, covarV = NULL){
  
  get.g=function(pair = NULL, phenotype, covarV){
    
    if(is.null(pair) || unique(pair) == "overall"){pair = NULL}
    
    if(is.null(pair)){
      pair.name <- "overall"
    }else{
      pair.name <- paste(pair, collapse = ",")
    }
    
    if(class(kin.obj) == "matrix"){
      K <- kin.obj
    }else{
      kin.mat.locale <- which(names(kin.obj) == pair.name)
      K <- kin.obj[[kin.mat.locale]]
    }
    
    #if we are correcting the covariate only don't put it in the model
    if(is.null(covarV) || is.null(pair)){
      model = regress::regress(as.vector(phenotype)~1,~K, pos = c(TRUE, TRUE))	
    }else{
      model = regress::regress(as.vector(phenotype)~covarV, ~K, pos = c(TRUE, TRUE))
    }
    
    #This err.cov is the same as err.cov in Dan's code using estVC
    #err.cov = summary(model)$sigma[1]*K+summary(model)$sigma[2]*diag(nrow(K))
    err.cov = model$sigma[1]*K+model$sigma[2]*diag(nrow(K))
    
    eW = eigen(err.cov, symmetric = TRUE)
    if(min(eW$values) < 0 && abs(min(eW$values)) > sqrt(.Machine$double.eps)){
    }else{
      eW$values[eW$values <= 0] = Inf
    } 
    err.cov = eW$vector %*% diag(eW$values^-0.5) %*% t(eW$vector)
    
    new.pheno <- err.cov %*% phenotype
    
    if(length(dim(geno)) == 3){
      l.geno <- lapply(1:dim(geno)[2], function(x) err.cov %*% geno[,x,]); #hist(new.geno)
      new.geno <- array(NA, dim = dim(geno))
      for(i in 1:length(l.geno)){
        new.geno[,i,] <- l.geno[[i]]
      }
      dimnames(new.geno) <- dimnames(geno)
    }else{
      new.geno <- err.cov %*% geno
      dimnames(new.geno) <- dimnames(geno)				
    }
    
    if(!is.null(covarV)){
      new.covar <- err.cov %*% covarV
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
