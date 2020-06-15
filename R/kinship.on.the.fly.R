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
      full.kin <- kin.obj
    }else{
      kin.mat.locale <- which(names(kin.obj) == pair.name)
      full.kin <- kin.obj[[kin.mat.locale]]
    }
    
    #also remove individuals with NAs in the phenotype
    #or covariates
    is.na.pheno <- which(is.na(phenotype))
    is.na.covar <- unique(which(is.na(covarV), arr.ind = TRUE)[,1])
    all.na <- unique(c(is.na.pheno, is.na.covar))
    not.na <- setdiff(1:length(phenotype), all.na)
    no.na.ind <- rownames(phenotype)[not.na]
    common.ind <- intersect(no.na.ind, colnames(full.kin))

    kin.locale <- match(common.ind, colnames(full.kin))  
    K <- full.kin[kin.locale,kin.locale]
    pheno.locale <- match(common.ind, rownames(phenotype))

    #if we are correcting the covariate only don't put it in the model
    if(is.null(covarV) || is.null(pair)){
      model = regress(as.vector(phenotype[pheno.locale])~1,~K, pos = c(TRUE, TRUE))	
    }else{
      model = regress(as.vector(phenotype)[pheno.locale]~covarV[pheno.locale,], ~K, pos = c(TRUE, TRUE))
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
