#' Perform linear regression on multi-allele markers.
#' 
#' This function performs the multi-allele version of 
#' linear regression. It is used in \code{\link{singlescan}}
#' and \code{\link{one_singlescanDO}}. It performs marker-by-marker
#' linear regressions for each trait and adjusts for covariates.
#' It collects parameters from the linear models and returns
#' for downstream use.
#' 
#' @param phenotype A phenotype vector
#' @param genotype genotypes of the marker being tested. 
#' If this is a vector, it will be converted into a one-column 
#' matrix. If it is a matrix, each column represents an allele,
#' and each row represents an individual.
#' @param covar_table The covariate matrix, with each covariate
#' in a column, and each individual in a row.
#' @param ph_family a character string for a description of the 
#' error distribution. Can be either "gaussian" or "binomial"
#' @param ref_col The column belonging to the reference allele.
#' The reference allele is removed from the genotype matrix so
#' that the matrix is linearly independent. There are only n-1
#' degrees of freedom in the genotype matrix, where n is the number
#' of alleles.
#' 
#' @return a list with "stats", "pval", "score"
#' stats is a matrix holding the t statistics and slopes (beta coefficients)
#' from the linear model.
#' pval holds the p value for the marker overall
#' and score holds the test statistic for the marker overall.
#'
#' @importFrom stats anova
#' 
#' @keywords internal
#' @export
#' 
get_stats_multiallele <- function(phenotype, genotype, covar_table, ph_family, ref_col){
  
  if(is.null(dim(genotype))){
    genotype <- matrix(genotype, ncol = 1)
  }
  
  #figure out if we are looking at a covariate
  covar_which <- check_geno(genotype, covar_table)
  
  #remove the column belonging to the reference allele
  if(ncol(genotype) > 1){
    geno_mat <- genotype[,-ref_col,drop=FALSE]
  }else{
    geno_mat <- genotype	
  }
  
  if(ncol(geno_mat) == 1){
    param_names <- "geno_mat"
  }else{
    param_names <- paste("geno_mat", colnames(geno_mat), sep = "")	
  }
  
  if(!is.null(covar_table)){
    covar_mat <- covar_table
  }else{
    covar_mat <- matrix(NA, nrow = length(phenotype), ncol = 0)
  }			
  
  if(length(covar_which) > 0){ #if we are looking at a covariate, take it out of the covariate matrix before testing
    corrected_covar <- covar_mat[,-covar_which, drop = FALSE]
  }else{
    corrected_covar <- covar_mat
  }
  
  #now do a check to see how many alleles of each
  #type are present
  bad_alleles <- NULL
  if(!is.null(dim(genotype)) && ncol(genotype) > 1){
    allele_counts <- apply(geno_mat, 2, function(x) length(unique(x)))
    bad_alleles <- as.vector(which(allele_counts < 2))
  }
  
  if(ncol(corrected_covar) == 0){
    model <- glm(phenotype~geno_mat, family = ph_family)
    model_coef <- coef(summary(model))
    coef_locale <- match(param_names, rownames(model_coef))
    
    slopes <- model_coef[coef_locale,1]
    ses <- model_coef[coef_locale,2]
    t_stats <- model_coef[coef_locale,3]
    locus_pval <- model_coef[coef_locale,4]
    
    if(ph_family == "gaussian"){
      locus_test <- anova(model, test = "F")
      locus_score <- as.matrix(locus_test)[2,4]
    }else{
      locus_test <- anova(model, test = "Chisq")
      locus_score <- as.matrix(locus_test)[2,4]
    }
    
  }else{
    model <- glm(phenotype~corrected_covar+geno_mat, family = ph_family)
    model_coef <- coef(summary(model))
    coef_locale <- match(param_names, rownames(model_coef))
    
    if(length(model_coef) < (ncol(corrected_covar)+ncol(geno_mat))){
      stop("I cannot fit coefficients to all parameters. Please check the covariate matrix.")
    }
    
    is_na <- which(is.na(cbind(corrected_covar, geno_mat)))
    if(length(is_na) > 0){
      not_na <- Reduce("intersect", apply(cbind(corrected_covar, geno_mat), 2, function(x) which(!is.na(x))))
    }else{
      not_na <- 1:nrow(geno_mat)	
    }
    full_model <- glm(phenotype[not_na]~cbind(corrected_covar, geno_mat)[not_na,], family = ph_family)
    bare_model <- glm(phenotype[not_na]~corrected_covar[not_na,], family = ph_family)
    if(ph_family == "gaussian"){
      comparison <- anova(bare_model, full_model, test = "F")
      locus_pval <- comparison$"Pr(>F)"[2]
      locus_score <- comparison$F[2]
    }else{
      comparison <- anova(bare_model, full_model, test = "Chisq")	
      locus_pval <- comparison$"Pr(>Chi)"[2]
      locus_score <- comparison$"Deviance"[2]
    }
    
    
    #if we are testing one of the covariates, the structure of the results is different
    #only fill in results for the A allele, all the others are the same
    if(length(covar_which) > 0){ 
      slopes <- c(rep(model_coef[coef_locale, 1], dim(geno_mat)[2]))
      ses <- c(rep(model_coef[coef_locale, 2], dim(geno_mat)[2]))
      t_stats <- c(rep(model_coef[coef_locale, 3], dim(geno_mat)[2]))
      p_vals <- c(rep(model_coef[coef_locale, 4], dim(geno_mat)[2]))
    }else{
      slopes <- as.vector(model_coef[coef_locale,1])
      ses <- as.vector(model_coef[coef_locale,2])
      t_stats <- as.vector(model_coef[coef_locale,3])
      p_vals <- as.vector(model_coef[coef_locale,4])
    } #end case for whether there are covariates specified
    
  } #end check for covariates
  
  
  #if there were some missing alleles, 0 out their entries in the vectors
  if(length(bad_alleles) > 0){
    # print(bad_alleles)
    # return(genotype)
    slopes[bad_alleles] <- 0
    ses[bad_alleles] <- 0
    t_stats[bad_alleles] <- 0
    # p_vals[bad_alleles] <- 1
  }
  #put together all the statistics we want to keep
  #we keep the slope, the t statistic, and the 
  results_table <- matrix(c(slopes, t_stats), byrow = TRUE, nrow = 2)
  colnames(results_table) <- colnames(geno_mat)
  rownames(results_table) <- c("slope", "t_stat")
  results <- list(results_table, locus_pval, locus_score); names(results) <- c("stats", "pval", "score")
  return(results)
}