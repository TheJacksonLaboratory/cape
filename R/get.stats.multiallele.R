#' This function performs regressions and
#' collects the relevant statistics.
#' Take out markers on the sex chromosomes
#' 
#' This function removes any markers in the geno.obj on the sex
#' chromosomes as well as invariant markers.
#'
#' @param phenotype kinship corrected phenotype data
#' @param genotype kinship corrected genotype data
#' @param covar.table kinship corrected covariate data
#' @param ph.family a character string for a description of the error distribution (see? "family" param in the glm function)
#' @param ref.col
#' 
#' @return a list with "stats", "pval", "score"
#'
#' @export
get.stats.multiallele <- function(phenotype, genotype, covar.table, ph.family, ref.col){
  
  if(is.null(dim(genotype))){
    genotype <- matrix(genotype, ncol = 1)
  }
  
  #figure out if we are looking at a covariate
  covar.which <- check.geno(genotype, covar.table)
  
  #remove the column belonging to the reference allele
  if(ncol(genotype) > 1){
    geno.mat <- genotype[,-ref.col,drop=FALSE]
  }else{
    geno.mat <- genotype	
  }
  
  if(ncol(geno.mat) == 1){
    param.names <- "geno.mat"
  }else{
    param.names <- paste("geno.mat", colnames(geno.mat), sep = "")	
  }
  
  if(!is.null(covar.table)){
    covar.mat <- covar.table
  }else{
    covar.mat <- matrix(NA, nrow = length(phenotype), ncol = 0)
  }			
  
  if(length(covar.which) > 0){ #if we are looking at a covariate, take it out of the covariate matrix before testing
    corrected.covar <- covar.mat[,-covar.which, drop = FALSE]
  }else{
    corrected.covar <- covar.mat
  }
  
  #now do a check to see how many alleles of each
  #type are present
  bad.alleles <- NULL
  if(!is.null(dim(genotype)) && ncol(genotype) > 1){
    allele.counts <- apply(geno.mat, 2, function(x) length(unique(x)))
    bad.alleles <- as.vector(which(allele.counts < 2))
  }
  
  if(ncol(corrected.covar) == 0){
    model <- glm(phenotype~geno.mat, family = ph.family)
    model.coef <- coef(summary(model))
    coef.locale <- match(param.names, rownames(model.coef))
    
    slopes <- model.coef[coef.locale,1]
    ses <- model.coef[coef.locale,2]
    t.stats <- model.coef[coef.locale,3]
    locus.pval <- model.coef[coef.locale,4]
    
    if(ph.family == "gaussian"){
      locus.test <- anova(model, test = "F")
      locus.score <- as.matrix(locus.test)[2,4]
    }else{
      locus.test <- anova(model, test = "Chisq")
      locus.score <- as.matrix(locus.test)[2,4]
    }
    
  }else{
    model <- glm(phenotype~corrected.covar+geno.mat, family = ph.family)
    model.coef <- coef(summary(model))
    coef.locale <- match(param.names, rownames(model.coef))
    
    if(length(model.coef) < (ncol(corrected.covar)+ncol(geno.mat))){
      stop("I cannot fit coefficients to all parameters. Please check the covariate matrix.")
    }
    
    is.na <- which(is.na(cbind(corrected.covar, geno.mat)))
    if(length(is.na) > 0){
      not.na <- Reduce("intersect", apply(cbind(corrected.covar, geno.mat), 2, function(x) which(!is.na(x))))
    }else{
      not.na <- 1:nrow(geno.mat)	
    }
    full.model <- glm(phenotype[not.na]~cbind(corrected.covar, geno.mat)[not.na,], family = ph.family)
    bare.model <- glm(phenotype[not.na]~corrected.covar[not.na,], family = ph.family)
    if(ph.family == "gaussian"){
      comparison <- anova(bare.model, full.model, test = "F")
      locus.pval <- comparison$"Pr(>F)"[2]
      locus.score <- comparison$F[2]
    }else{
      comparison <- anova(bare.model, full.model, test = "Chisq")	
      locus.pval <- comparison$"Pr(>Chi)"[2]
      locus.score <- comparison$"Deviance"[2]
    }
    
    
    #if we are testing one of the covariates, the structure of the results is different
    #only fill in results for the A allele, all the others are the same
    if(length(covar.which) > 0){ 
      slopes <- c(rep(model.coef[coef.locale, 1], dim(geno.mat)[2]))
      ses <- c(rep(model.coef[coef.locale, 2], dim(geno.mat)[2]))
      t.stats <- c(rep(model.coef[coef.locale, 3], dim(geno.mat)[2]))
      p.vals <- c(rep(model.coef[coef.locale, 4], dim(geno.mat)[2]))
    }else{
      slopes <- as.vector(model.coef[coef.locale,1])
      ses <- as.vector(model.coef[coef.locale,2])
      t.stats <- as.vector(model.coef[coef.locale,3])
      p.vals <- as.vector(model.coef[coef.locale,4])
    } #end case for whether there are covariates specified
    
  } #end check for covariates
  
  
  #if there were some missing alleles, 0 out their entries in the vectors
  if(length(bad.alleles) > 0){
    # print(bad.alleles)
    # return(genotype)
    slopes[bad.alleles] <- 0
    ses[bad.alleles] <- 0
    t.stats[bad.alleles] <- 0
    # p.vals[bad.alleles] <- 1
  }
  #put together all the statistics we want to keep
  #we keep the slope, the t statistic, and the 
  results.table <- matrix(c(slopes, t.stats), byrow = TRUE, nrow = 2)
  colnames(results.table) <- colnames(geno.mat)
  rownames(results.table) <- c("slope", "t.stat")
  results <- list(results.table, locus.pval, locus.score); names(results) <- c("stats", "pval", "score")
  return(results)
}