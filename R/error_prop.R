#' Estimate Errors of Regression Coefficients
#' 
#' This function uses error propagation formulas for 
#' quantities computed from regression coefficients to 
#' estimate the error for all regression coefficients.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param pairscan_obj a pairscan object from \code{\link{pairscan}}
#' @param perm A logical value to indicate whether error propagation should be 
#' performed on the test statistics (FALSE) or the permuted test statistics (TRUE).
#' @param verbose A logical value to indicate whether the progress of the function 
#' should be printed to the screen.
#' @param run_parallel boolean, default = FALSE
#' @param n_cores The number of cores to use if run_parallel is TRUE, default = 4
#' @param just_m If TRUE only the m12 and m21 values are calculated. If FALSE, the
#' default, the standard deviations are also calculated.
#' 
#' @return This function returns the data object with a new list element: var_to_var_influences 
#' if perm is set to FALSE and var_to_var_influences_perm if perm is set to TRUE. These tables 
#' include the errors calculated for the marker1 to marker2 (m21) influences as well as the
#' marker2 to marker1 (m12) influences. These results are used by \code{\link{calc_p}} to
#' calculate empirical p values.
#' 
#' @import parallel
#' @import foreach
#' @importFrom doParallel registerDoParallel 
#' 
#' @examples 
#' \dontrun{
#' #run error propagateion on test statistics and 
#' #permuted test statistics
#' data_obj <- error_prop(data_obj, pairscan_obj, perm = TRUE)
#' data_obj <- error_prop(data_obj, pairscan_obj, perm = FALSE)
#' }
#' @export
error_prop <- function (data_obj, pairscan_obj, perm = FALSE, verbose = FALSE,
                        run_parallel = FALSE, n_cores = 4, just_m = FALSE) {
  
  if(!run_parallel){n_cores = 1}
  
  scan_two_results = NULL
  
  if(perm){
    scan_two_results <- pairscan_obj$pairscan_perm
  }else{
    scan_two_results <- pairscan_obj$pairscan_results
  }
  
  p = NULL #for appeasing R CMD check
  
  results_obj <- vector(mode = "list", length = 2)
  names(results_obj) <- c("var_to_var_influences_perm")
  
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
  get_beta_mat <- function(marker_pair_number, scan_two_results){
    #the beta matrix is composed of the coefficients from each pairwise
    #marker model (except for the interaction coefficient)
    #we use all markers in each row and set the non-covariate entries to 0
    num_pheno <- length(scan_two_results)
    beta_mat <- sapply(scan_two_results, function(x) as.numeric(x[[1]][marker_pair_number,(dim(scan_two_results[[1]][[1]])[2]-2):(dim(scan_two_results[[1]][[1]])[2])]))
    rownames(beta_mat) <- c("marker1", "marker2", "interaction")
    return(beta_mat)	
  }
  
  
  get_se_mat <- function(marker_pair_number, scan_two_results){
    #the beta matrix is composed of the coefficients from each pairwise
    #marker model (except for the interaction coefficient)
    #we use all markers in each row and set the non-covariate entries to 0
    num_pheno <- length(scan_two_results)
    se_mat <- sapply(scan_two_results, function(x) as.numeric(x[[2]][marker_pair_number,(dim(scan_two_results[[1]][[1]])[2]-2):(dim(scan_two_results[[1]][[1]])[2])]))
    rownames(se_mat) <- c("marker1", "marker2", "interaction")
    return(se_mat)	
  }
  
  get_cov_mat <- function(marker_pair_number, num_pheno){
    #the variance-covariance matrix is block diagonal, and 
    #contains three times the number of rows and columns
    #as scanned traits. The blocks contain the variance-
    #covariance matrix from each trait for the marker pair
    cov_mat <- matrix(0, num_pheno*3, num_pheno*3)
    pheno_num <- 1:num_pheno
    start_row_col <- (pheno_num*3)-2
    end_row_col <- (pheno_num*3)
    for(ph in 1:num_pheno){
      cov_mat[start_row_col[ph]:end_row_col[ph], start_row_col[ph]:end_row_col[ph]] <- matrix(as.numeric(scan_two_results[[ph]][[3]][marker_pair_number,]), ncol = 3)
    }	
    return(cov_mat)
  }
  
  #check to see if phenotypes or eigentraits were scanned
  pheno_names <- names(pairscan_obj$pairscan_results)	
  pheno_check <- match(pheno_names, colnames(data_obj$pheno))
  if(length(which(!is.na(pheno_check))) == 0){ #if we scanned eigentraits
    num_pheno <- dim(data_obj$ET)[2] #the number of phenotypes
    names_pheno <- colnames(data_obj$ET)
  }else{
    num_pheno <- dim(data_obj$pheno)[2] #the number of phenotypes
    names_pheno <- colnames(data_obj$pheno)
  }				
  
  get_pair_coeffs <- function(marker_pair_number, scan_two_results, marker_mat){
    beta_main <- t(get_beta_mat(marker_pair_number, scan_two_results)) ### Extract Main effect and interactions
    non_zero <- which(beta_main[1:2,] != 0)
    na_locale <- which(is.na(beta_main))
    if(length(non_zero) > 0 && length(na_locale) == 0){
      num_pheno <- length(scan_two_results)
      beta_se <- t(get_se_mat(marker_pair_number, scan_two_results)) ### Extract Main effect and interactions
      beta_cov <- get_cov_mat(marker_pair_number, num_pheno) ### Extract Covars
      if(!just_m){
        inf_coeffs <- calc_delta_errors(markers = marker_mat[marker_pair_number,], 
        beta_m = beta_main, se = beta_se, beta_cov)
      }else{
        inf_coeffs <- calc_m(markers = marker_mat[marker_pair_number,], beta_m = beta_main)	
      }
    }else{
      inf_coeffs <- NULL	
    }
    return(inf_coeffs)
  }
  
  get_multi_pair_coeffs <- function(marker_pair_number_v, scan_two_results, marker_mat){
    result <- do.call("rbind", lapply(marker_pair_number_v, function(x) get_pair_coeffs(x, scan_two_results, marker_mat)))
    return(result)
  }
  
  #====================================================================================
  #end internal functions
  #====================================================================================
  
  ### For all marker pairs calculate activity and IC
  
  if(perm){
    if(is.null(pairscan_obj$pairscan_perm)){
      stop("pairscan() with permutations must be run before error_prop()")
    }
    marker_mat <-  pairscan_obj$pairscan_perm[[1]][[1]][,1:2]#get all the pairs that were tested in the pair scan
    scan_two_results <- pairscan_obj$pairscan_perm  #the coefficient matrices from the 2D scan			
  }else{
    if(is.null(pairscan_obj$pairscan_results)){
      stop("pairscan() must be run before error_prop()")
    }
    marker_mat <- pairscan_obj$pairscan_results[[1]][[1]][,1:2] #a matrix listing the names of used marker combinations
    scan_two_results <- pairscan_obj$pairscan_results  #the coefficient matrices from the 2D scan
  }
  
  
  colnames(marker_mat) <- c("marker1", "marker2")
  n_pairs <- length(marker_mat[,1]) #number of pairs of genes
  
  chunked_pairs <- chunkV(1:n_pairs, n_cores)
  
  if (run_parallel) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    cape_dir_full <- find.package("cape")
    cape_dir <- str_replace(cape_dir_full,"cape_pkg/cape","cape_pkg")
    clusterExport(cl, "cape_dir", envir=environment())
    clusterEvalQ(cl, .libPaths(cape_dir))
    # , .packages = 'cape'
    influence_coeffs <- foreach(p = c(1:length(chunked_pairs)), .combine = "rbind", .export = c("calc_delta_errors", "pseudoinverse", "propagate", "calc_m")) %dopar% {
      get_multi_pair_coeffs(chunked_pairs[[p]], scan_two_results, marker_mat)
    }				
    stopCluster(cl)
    
  } else {
    
    influence_coeffs <- c()
    index <- 1:length(chunked_pairs)
    for (p in index) {
      influence_coeffs <- rbind(influence_coeffs, get_multi_pair_coeffs(chunked_pairs[[p]], scan_two_results, marker_mat))
    }
    
  }
  
  if(!just_m){
    colnames(influence_coeffs) <- c("marker1","marker2","m12","m12_std_dev","m21","m21_std_dev")
  }else{
    colnames(influence_coeffs) <- c("marker1","marker2","m12","m21")	
  }
  
  if(perm){
    data_obj$var_to_var_influences_perm <- influence_coeffs
  }else{
    data_obj$var_to_var_influences <- influence_coeffs			
  }		
  
  if(verbose){
    cat("\n") #make sure the prompt is on the next line at the end of everything
  }
  
  return(data_obj)
}
