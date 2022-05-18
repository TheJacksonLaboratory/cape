#' This is an internal function to run a single pairscan
#' It is used both to do the actual pairscan 
#' (\code{\link{pairscan_kin}} and \code{\link{pairscan_noKin}}), 
#' as well as to do the permutations of the pairscan
#' \code{\link{pairscan_null}}).
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param phenotype_vector A vector of trait values
#' @param genotype_matrix A matrix of genotypes for markers to
#' be tested
#' @param int the interaction term added to the linear model
#' when the kinship correction is being used. This term is
#' calculated in \code{\link{pairscan_kin}}.
#' @param covar_vector a vector or matrix of covariates to
#' be used.
#' @param paired_markers a two-column matrix indicating which
#' marker pairs should be tested. The pairs are assigned in
#' \code{\link{pairscan}} by \code{\link{get_pairs_for_pairscan}}.
#' They are checked for pairwise correlations before being sent
#' to the pairscan.
#' @param n_perm the number of permutations to be performed.
#' @param run_parallel a logical value indicating whether to 
#' use parallel computing
#' @param verbose a logical value indicating whether progress should
#' be printed to the screen.
#' @param n_cores the number of CPUs to use if run_parallel is TRUE.
#'
#' @return This function returns a list with two slots: 
#' pairscan_results and pairscan_perm
#' Each of these elements is also a list containing effect
#' sizes, standard errors, and covariance matrices for the
#' pairwise tests.
#' 
#' @import parallel
#' @import foreach
#' @importFrom Matrix rankMatrix
#' @importFrom doParallel registerDoParallel
#' @importFrom stats coefficients vcov
#' @keywords internal
#' 
one_pairscan_parallel <- function(data_obj, phenotype_vector, genotype_matrix, 
  int = NULL, covar_vector = NULL, paired_markers, n_perm = 0, run_parallel = FALSE, 
  verbose = FALSE, n_cores = 4){
  
  if(!run_parallel){n_cores = 1}		
  
  m = p = NULL #for appeasing R CMD check
  covar_names <- get_covar(data_obj)$covar_names  # don't change this to underscore notation!
  
  #============================================================================
  # check to see that the covariates are not redundant and are linearly independent
  #============================================================================
  use_covars <- as.logical(length(covar_vector) > 0)
  if(use_covars > 0){
    cov_mat <- covar_vector
    
    #remove the NAs and check the matrix for rank
    not_na_locale <- which(!is.na(rowSums(cov_mat)))
    no_na_cov <- cov_mat[not_na_locale,,drop=FALSE]
    
    design_cov <- cbind(rep(1, dim(no_na_cov)[1]), no_na_cov)
    rank_cov <- rankMatrix(design_cov)
    if(rank_cov[[1]] < dim(design_cov)[2]){
      stop("The covariate matrix does not appear to be linearly independent.\nIf you are using dummy variables for groups, leave one of the groups out.")
    }
    
    cor_mat <- cor(no_na_cov)
    diag(cor_mat) <- 0
    perfect_cor <- which(abs(signif(cor_mat, 2)) == 1)
    if(length(perfect_cor) > 0){
      stop("Check the covariates. There appears to be at least one pair of redundant covariates.")
    }
    num_covar <- dim(covar_vector)[2]	
  }else{
    num_covar = 0
  }			
  #============================================================================
  
  
  #============================================================================
  #internal functions
  #============================================================================
  
  
  get_model_results <- function(marker_names, m1, m2, int_term = NULL, testing_covar = FALSE){
    #if we are testing a covariate, pull it out of the 
    #covariate.table
    if(testing_covar){
      covar_locale <- which(covar_names %in% marker_names)
      if(length(covar_locale) > 0){
        new_covar_vector <- covar_vector[,-covar_locale,drop=FALSE]
      }else{
        new_covar_vector <- covar_vector	
      }
    }else{
      new_covar_vector <- covar_vector
    }
    
    if(is.null(new_covar_vector) || dim(new_covar_vector)[2] == 0){
      new_covar_vector <- NULL
    }
    
    design_mat <- cbind(rep(1, length(m1)), new_covar_vector, m1, m2)
    
    missing_rows <- which(is.na(rowSums(design_mat)))
    if(length(missing_rows) > 0){
      design_mat <- design_mat[-missing_rows,]
    }
    rank_check <- rankMatrix(design_mat)[[1]]
    
    if(rank_check < dim(design_mat)[2]){
      return(NULL)
    }
    
    #Do the linear regression with the covariates, the two markers 
    #individually, and the interaction between the two markers
    #put the covars first so the marker effects come last
    if(!is.null(new_covar_vector)){
      if(is.null(int_term)){
        model <- lm(phenotype_vector ~ new_covar_vector + m1 + m2 + m1:m2, model = FALSE, x = FALSE, y = FALSE, qr = TRUE)
      }else{
        model <- lm(phenotype_vector ~ 0 + new_covar_vector + m1 + m2 + int_term, model = FALSE, x = FALSE, y = FALSE, qr = TRUE)	
      }
    }else{
      if(is.null(int_term)){
        model <- lm(phenotype_vector ~ m1 + m2 + m1:m2, model = FALSE, x = FALSE, y = FALSE, qr = TRUE)	
      }else{
        model <- lm(phenotype_vector ~ 0 + m1 + m2 + int_term, model = FALSE, x = FALSE, y = FALSE, qr = TRUE)		
      }
    }
    model_summ <- summary(model)
    
    #occasionally permutations result in non linearly dependent matrices.
    #if this is the case, return NULL. This triggers the permutation
    #script to generate another permutation.
    if(length(which(is.na(coefficients(model)))) > 0){ 
      return(NULL)
    }
    
    #take the last 3 terms
    keep_results <- tail(1:length(coef(model)), 3)
    model_effects <- model_summ$coefficients[keep_results,"Estimate"]
    model_se <- model_summ$coefficients[keep_results,"Std. Error"]
    
    #calculate the covariance matrix of the model parameters
    model_cov_mat <- vcov(model)
    dim_mat <- dim(model_cov_mat)[1] #get the dimensions of the matrix. We want the last three rows and last three columns (the covariance matrix for m1, m2, and m1:m2)
    cov_mat <- model_cov_mat[(dim_mat-2):dim_mat, (dim_mat-2):dim_mat]
    model_cov_results <- as.vector(cov_mat)
    
    results <- list(model_effects, model_se, model_cov_results)
    names(results) <- c("model_effects", "model_se", "model.cov")
    return(results)
  }
  
  get_pair_results <- function(m_pair){
    #check the markers for a covariate
    covar_check <- which(m_pair %in% covar_names)
    if(length(covar_check) > 0){
      testing_covar = TRUE
    }else{
      testing_covar = FALSE	
    }
    
    #get the marker identities
    marker1 <- genotype_matrix[,as.character(m_pair[1])]
    marker2 <- genotype_matrix[,as.character(m_pair[2])]
    
    if(is.null(int)){
      marker_pair_results <- get_model_results(marker_names = m_pair, m1 = marker1, m2 = marker2, testing_covar = testing_covar)
    }else{
      marker_pair_results <- get_model_results(marker_names = m_pair, m1 = marker1, m2 = marker2, int_term = int, testing_covar = testing_covar)
    }
    
    return(marker_pair_results)				
  }
  
  
  get_multi_pair_results <- function(m_pair_v){
    result <- lapply(m_pair_v[[1]], function(x) get_pair_results(m_pair = paired_markers[x,]))
    return(result)
  }
  
  
  
  one_perm <- function(perm.num){
    m_pair <- paired_markers[sample(1:dim(paired_markers)[1], 1),]
    rnd_order <- sample(1:dim(genotype_matrix)[1], dim(genotype_matrix)[1])
    marker1 <- genotype_matrix[rnd_order,as.character(m_pair[1])]
    marker2 <- genotype_matrix[rnd_order,as.character(m_pair[2])]
    
    if(is.null(int)){
      marker_pair_results <- get_model_results(marker_names = m_pair, m1 = marker1, m2 = marker2)
    }else{
      marker_pair_results <- get_model_results(marker_names = m_pair, m1 = marker1, m2 = marker2, int_term = int)	
    }
    
    marker_pair_results$"pair.used" <- m_pair
    return(marker_pair_results)				
  }
  #============================================================================
  #end of internal functions
  #============================================================================
  
  # each.iter.time <- rep(NA, n_pairs)
  #we can calculate the number of genes
  #from the input data
  # TODO remove the next line? n_pairs is not used in this script
  # n_pairs <- dim(paired_markers)[1]
  
  if(nrow(paired_markers) > n_cores){
    chunked_pairs <- chunkV(1:nrow(paired_markers), n_cores)
  }else{
    chunked_pairs <- list(1:nrow(paired_markers))
  }
  
  if (run_parallel) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    cape_dir_full <- find.package("cape")
    cape_dir <- gsub("cape_pkg/cape","cape_pkg", cape_dir_full)
    clusterExport(cl, "cape_dir", envir=environment())
    clusterEvalQ(cl, .libPaths(cape_dir))
    pair_results_list <- foreach(m = 1:length(chunked_pairs), .packages = 'cape', .export = c("phenotype_vector", "rankMatrix")) %dopar% {
      get_multi_pair_results(m_pair_v = chunked_pairs[m])
    }
    stopCluster(cl)
    
  } else {
    
    pair_results_list <- vector(mode = "list", length = length(chunked_pairs))
    for(ind in length(chunked_pairs)) {
      pair_results_list[[ind]] <- get_multi_pair_results(m_pair_v = chunked_pairs[ind])
    }
    
  }
  
  pair_results <- unlist(pair_results_list, recursive = FALSE)
  
  #Filter out the results with null values
  good_results <- which(unlist(lapply(pair_results, function(x) length(x$model_effects))) > 0)
  if(length(good_results) > 0){
    pair_results <- pair_results[good_results]
    all_model_effects <- matrix(unlist(lapply(pair_results, function(x) x$model_effects)), nrow = length(pair_results), byrow = TRUE)
    all_model_se <- matrix(unlist(lapply(pair_results, function(x) x$model_se)), nrow = length(pair_results), byrow = TRUE)
    all_model_cov <- matrix(unlist(lapply(pair_results, function(x) x$model.cov)), nrow = length(pair_results), byrow = TRUE)
    
    #assign column names to the results tables
    #the column names represent the possible beta
    #coefficients we can get. The intercept, all
    #possible covariates, marker1 and marker2,
    #and the interaction marker1:marker2
    column_names <- c("marker1", "marker2", "marker1:marker2")	
    colnames(all_model_effects) <- colnames(all_model_se) <- column_names
    
    # add the marker pair names to the results tables and name the columns
    marker_labels <- paired_markers[good_results,,drop=FALSE]
    colnames(marker_labels) <- c("marker_name1", "marker_name2")
    final_effects_table <- cbind(paired_markers[good_results,,drop=FALSE], all_model_effects)
    final_se_table <- cbind(paired_markers[good_results,,drop=FALSE], all_model_se)
    final_cov_table <- all_model_cov
    
    
    phenotype_results <- list(final_effects_table, final_se_table, final_cov_table)
    names(phenotype_results) <- c("pairscan_effects", "pairscan_se", "model_covariance")
  }else{
    phenotype_results <- NULL	
  }
  
  
  if(n_perm > 0){
    if(verbose){cat("\tCalculating permutations...\n")}
    
    if (run_parallel) {
      
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      cape_dir_full <- find.package("cape")
      cape_dir <- gsub("cape_pkg/cape","cape_pkg", cape_dir_full)
      clusterExport(cl, "cape_dir", envir=environment())
      clusterEvalQ(cl, .libPaths(cape_dir))
      perm_results <- foreach(p = 1:n_perm, .packages = 'cape', .export = c("phenotype_vector", "rankMatrix")) %dopar% {
        one_perm(p)
      }
      stopCluster(cl)
      
    } else {
      
      perm_results <- c()
      index <- 1:n_perm
      for (p in index) {
        perm_results <- rbind(perm_results, one_perm(p))
      }
    }
    
    #also make variables to hold the permutation results
    good_results_perm <- which(unlist(lapply(perm_results, function(x) length(x$model_effects))) > 0)
    perm_results <- perm_results[good_results_perm]			
    all_model_effects_perm <- matrix(unlist(lapply(perm_results, function(x) x$model_effects)), nrow = length(good_results_perm), byrow = TRUE)
    all_model_se_perm <- matrix(unlist(lapply(perm_results, function(x) x$model_se)), nrow = length(good_results_perm), byrow = TRUE)
    all_model_cov_perm <- matrix(unlist(lapply(perm_results, function(x) x$model.cov)), nrow = length(good_results_perm), byrow = TRUE)
    marker_pairs_used_perm <- matrix(unlist(lapply(perm_results, function(x) x$pair.used)), nrow = length(good_results_perm), byrow = TRUE)
    colnames(all_model_effects_perm) <- colnames(all_model_se_perm) <- column_names			
    rm(perm_results)
    
    colnames(marker_pairs_used_perm) <- c("marker1", "marker2")
    final_effects_table_perm <- cbind(marker_pairs_used_perm, all_model_effects_perm)
    final_se_table_perm <- cbind(marker_pairs_used_perm, all_model_se_perm)	
    
    final_cov_table_perm <- all_model_cov_perm
    phenotype_perm_results <- list(final_effects_table_perm, final_se_table_perm, final_cov_table_perm)
    names(phenotype_perm_results) <- c("pairscan_effects.perm", "pairscan_se.perm", "model_covariance.perm")
  }else{
    phenotype_perm_results <- NULL
  }
  
  final_results <- list(phenotype_results, phenotype_perm_results)
  names(final_results) <- c("pairscan_results", "pairscan_perm")
  
  #if(verbose){cat("\n")} #make sure the prompt is on the next line at the end of everything
  
  return(final_results)
  
}
