#' Estimate Errors of Regression Coefficients
#' 
#' This function uses error propagation formulas for quantities computed from 
#' regression coefficients to estimate the error for all regression coefficients.
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param pairscan.obj a pariscan object from \code{\link{pairscan}}
#' @param perm A logical value to indicate whether error propagation should be 
#' performed on the test statistics (FALSE) or the permuted test statistics (TRUE).
#' @param verbose A logical value to indicate whether the progress of the function 
#' should be printed to the screen.
#' @param run.parallel boolean, default = TRUE
#' @param n.cores integer, default = 4
#' @param just.m 
#' 
#' @return This function returns the data object with a new list element: var_to_var_influences 
#' if perm is set to FALSE and var.to.var.influences.perm if perm is set to TRUE. These tables 
#' include the errors calculated for the marker1 to marker2 influences as well as the marker2 
#' to marker1 influences. These results are used by \code{\link{calc.p}} to calculate empirical 
#' p values.
#' 
error.prop <- function (data.obj, pairscan.obj, perm = FALSE, verbose = FALSE,
                        run.parallel = FALSE, n.cores = 4, just.m = FALSE) {
  
  if(!run.parallel){n.cores = 1}
  
  scan.two.results = NULL
  
  if(perm){
    scan.two.results <- pairscan.obj$pairscan.perm
  }else{
    scan.two.results <- pairscan.obj$pairscan.results
  }
  
  p = NULL #for appeasing R CMD check
  
  results.obj <- vector(mode = "list", length = 2)
  names(results.obj) <- c("var.to.var.influences.perm")
  
  if(verbose){
    if(perm){
      cat("\nCalculating error propagation of permuted coefficients.\n")
    }else{
      cat("\nCalculating error propagation of coefficients.\n")	
    }
  }
  
  #====================================================================================
  #begin internal functions
  #====================================================================================
  get.beta.mat <- function(marker.pair.number, scan.two.results){
    #the beta matrix is composed of the coefficients from each pairwise
    #marker model (except for the interaction coefficient)
    #we use all markers in each row and set the non-covariate entries to 0
    num.pheno <- length(scan.two.results)
    beta.mat <- sapply(scan.two.results, function(x) as.numeric(x[[1]][marker.pair.number,(dim(scan.two.results[[1]][[1]])[2]-2):(dim(scan.two.results[[1]][[1]])[2])]))
    rownames(beta.mat) <- c("marker1", "marker2", "interaction")
    return(beta.mat)	
  }
  
  
  get.se.mat <- function(marker.pair.number, scan.two.results){
    #the beta matrix is composed of the coefficients from each pairwise
    #marker model (except for the interaction coefficient)
    #we use all markers in each row and set the non-covariate entries to 0
    num.pheno <- length(scan.two.results)
    se.mat <- sapply(scan.two.results, function(x) as.numeric(x[[2]][marker.pair.number,(dim(scan.two.results[[1]][[1]])[2]-2):(dim(scan.two.results[[1]][[1]])[2])]))
    rownames(se.mat) <- c("marker1", "marker2", "interaction")
    return(se.mat)	
  }
  
  get.cov.mat <- function(marker.pair.number, num.pheno){
    #the variance-covariance matrix is block diagonal, and 
    #contains three times the number of rows and columns
    #as scanned traits. The blocks contain the variance-
    #covariance matrix from each trait for the marker pair
    cov.mat <- matrix(0, num.pheno*3, num.pheno*3)
    pheno.num <- 1:num.pheno
    start.row.col <- (pheno.num*3)-2
    end.row.col <- (pheno.num*3)
    for(ph in 1:num.pheno){
      cov.mat[start.row.col[ph]:end.row.col[ph], start.row.col[ph]:end.row.col[ph]] <- matrix(scan.two.results[[ph]][[3]][marker.pair.number,], ncol = 3)
    }	
    return(cov.mat)
  }
  
  #check to see if phenotypes or eigentraits were scanned
  pheno.names <- names(pairscan.obj$pairscan.results)	
  pheno.check <- match(pheno.names, colnames(data.obj$pheno))
  if(length(which(!is.na(pheno.check))) == 0){ #if we scanned eigentraits
    num.pheno <- dim(data.obj$ET)[2] #the number of phenotypes
    names.pheno <- colnames(data.obj$ET)
  }else{
    num.pheno <- dim(data.obj$pheno)[2] #the number of phenotypes
    names.pheno <- colnames(data.obj$pheno)
  }				
  
  get.pair.coeffs <- function(marker.pair.number, scan.two.results, marker.mat){
    beta.main <- t(get.beta.mat(marker.pair.number, scan.two.results)) ### Extract Main effect and interactions
    non.zero <- which(beta.main[1:2,] != 0)
    na.locale <- which(is.na(beta.main))
    if(length(non.zero) > 0 && length(na.locale) == 0){
      num.pheno <- length(scan.two.results)
      beta.se <- t(get.se.mat(marker.pair.number, scan.two.results)) ### Extract Main effect and interactions
      beta.cov <- get.cov.mat(marker.pair.number, num.pheno) ### Extract Covars
      if(!just.m){
        inf.coeffs <- calc.delta.errors(markers = marker.mat[marker.pair.number,], beta.m = beta.main, se = beta.se, beta.cov)
      }else{
        inf.coeffs <- calc.m(markers = marker.mat[marker.pair.number,], beta.m = beta.main, beta.cov)	
      }
    }else{
      inf.coeffs <- NULL	
    }
    return(inf.coeffs)
  }
  
  get.multi.pair.coeffs <- function(marker.pair.number.v, scan.two.results, marker.mat){
    result <- do.call("rbind", lapply(marker.pair.number.v, function(x) get.pair.coeffs(x, scan.two.results, marker.mat)))
    return(result)
  }
  
  #====================================================================================
  #end internal functions
  #====================================================================================
  
  ### For all marker pairs calculate activity and IC
  
  if(perm){
    if(is.null(pairscan.obj$pairscan.perm)){
      stop("pairscan() with permutations must be run before error.prop()")
    }
    marker.mat <-  pairscan.obj$pairscan.perm[[1]][[1]][,1:2]#get all the pairs that were tested in the pair scan
    scan.two.results <- pairscan.obj$pairscan.perm  #the coefficient matrices from the 2D scan			
  }else{
    if(is.null(pairscan.obj$pairscan.results)){
      stop("pairscan() must be run before error.prop()")
    }
    marker.mat <- pairscan.obj$pairscan.results[[1]][[1]][,1:2] #a matrix listing the names of used marker combinations
    scan.two.results <- pairscan.obj$pairscan.results  #the coefficient matrices from the 2D scan
  }
  
  
  colnames(marker.mat) <- c("marker1", "marker2")
  n.pairs <- length(marker.mat[,1]) #number of pairs of genes
  
  chunked.pairs <- chunkV(1:n.pairs, n.cores)
  
  if (run.parallel) {
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    cape.dir <- paste(find.package("cape"),"/cape_pkg",sep="")
    parallel::clusterCall(cl, function() {.libPaths(cape.dir)})
    # , .packages = 'cape'
    influence.coeffs <- foreach::foreach(p = c(1:length(chunked.pairs)), .combine = "rbind", .export = c("calc.delta.errors", "pseudoinverse", "propagate", "calc.m")) %dopar% {
      get.multi.pair.coeffs(chunked.pairs[[p]], scan.two.results, marker.mat)
    }				
    parallel::stopCluster(cl)
    
  } else {
    
    influence.coeffs <- c()
    index <- 1:length(chunked.pairs)
    for (p in index) {
      influence.coeffs <- rbind(influence.coeffs, get.multi.pair.coeffs(chunked.pairs[[p]], scan.two.results, marker.mat))
    }
    
  }
  
  if(!just.m){
    colnames(influence.coeffs) <- c("marker1","marker2","m12","m12.std.dev","m21","m21.std.dev")
  }else{
    colnames(influence.coeffs) <- c("marker1","marker2","m12","m21")	
  }
  
  if(perm){
    data.obj$var_to_var_influences_perm <- influence.coeffs
  }else{
    data.obj$var_to_var_influences <- influence.coeffs			
  }		
  
  if(verbose){
    cat("\n") #make sure the prompt is on the next line at the end of everything
  }
  
  return(data.obj)
}
