#' Select marker pairs for pairscan
#' 
#' This function selects which marker pairs can be tested in the pair scan.
#' Even if all markers are linearly independent, some marker pairs may have
#' insufficient recombination between them to populate all genotype
#' combinations. Marker pairs for which genotype combinations have insufficient
#' numbers of individuals are not tested. This function determines which marker
#' pairs have sufficient representation in all genotype combinations. 
#'  
#' @param gene A two dimensional genotype matrix with rows containing 
#'   individuals and columns containing markers. Each entry is a value between
#'   0 and 1 indicating the genotype of each individual at each marker. 
#' @param covar_names A character vector indicating which covariates should
#' be tested.
#' @param min_per_genotype The minimum number of individuals allowable per 
#'   genotype. If for a given marker pair, one of the genotypes is 
#'   underrepresented, the marker pair is not tested. If this value is NULL,
#'   max_pair_cor must have a numeric value.
#' @param max_pair_cor A numeric value between 0 and 1 indicating the maximum
#'   Pearson correlation that two markers are allowed. If the correlation
#'   between a pair of markers exceeds this threshold, the pair is not tested.
#'   If this value is set to NULL, min_per_genotype must have a numeric value.
#' @param run_parallel A logical value indicating whether multiple 
#'   processors should be used.
#' @param n_cores The number of cores to be used if run_parallel is TRUE
#' @param verbose A logical value. If TRUE, the script prints a message to the
#'   screen to indicate that it is running. If FALSE, no message is printed.
#' 
#' @return This function returns a two-column matrix of marker pairs. This
#'   matrix is then used as an argument in \code{\link{one_pairscan_parallel}}, 
#'   \code{\link{pairscan_null_kin}}, \code{\link{pairscan_null}} and 
#'   \code{\link{pairscan}} to specify which marker pairs should be tested.
#' 
#' @details One and only one of min_per_genotype or max_pair_cor should be specified.
#' We recommend that if you have continuous genotype probabilities, you use max_pair_cor.
#' If both values are specified, this function will preferentially use max_pair_cor.
#' 
#' @import parallel
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' 
#' @export
get_pairs_for_pairscan <- function(gene, covar_names = NULL, max_pair_cor = NULL, 
	min_per_genotype = NULL, run_parallel = FALSE, n_cores = 4, verbose = FALSE){
  
  if(!run_parallel){n_cores = 1}
  
  p = NULL #for appeasing R CMD check
  
  if(is.null(max_pair_cor) && is.null(min_per_genotype)){
    stop("One of max_pair_cor or min_per_genotype should be set.")
  }
  
  if(!is.null(max_pair_cor) && !is.null(min_per_genotype)){
    warning("Only one of max_pair_cor or min_per_genotype should be set. Choosing only max_pair_cor.")
    min_per_genotype = NULL
  }
  
  if(verbose){
    cat("\nChecking marker pairs for genotype representation...\n")
  }
  
  # covariates should be paired with every marker
  # so only check pairs for genetic markers, and make
  # a separate list of pairs for covariates with
  # markers and each other
  all_covar_pairs <- NULL
  if(!is.null(covar_names)){
    covar_locale <- match(covar_names, colnames(gene))
    gene <- gene[,-covar_locale,drop=FALSE]
    covar_pairs <- pair_matrix(covar_names)
    covar_geno_pairs <- cbind(rep(colnames(gene), length(covar_names)), rep(covar_names, each = ncol(gene)))
    all_covar_pairs <- rbind(covar_geno_pairs, covar_pairs)
  }
  
  if(!is.null(max_pair_cor)){
    thresh_param <- max_pair_cor
    check_linkage <- function(m1,m2,thresh_param){
      pair_cor <- try(cor(m1, m2, use = "complete"), silent = TRUE)
      class_comp <- class(pair_cor)
      if(class_comp == "try-error" || pair_cor > max_pair_cor || is.na(pair_cor)) {
        return(FALSE) #pair failed check
      }else{
        return(TRUE) #pair passed check
      }
    }
  }
  
  
  if(!is.null(min_per_genotype)){
    thresh_param <- min_per_genotype
    check_linkage <- function(m1,m2,thresh_param){
      geno_table <- cbind.data.frame(as.factor(m1),as.factor(m2))
      colnames(geno_table) <- c("m1","m2")
      reps <- table(geno_table$m1,geno_table$m2)
      too_few <- which(reps < thresh_param)
      if(length(too_few) >= 1) {
        return(FALSE) #pair failed check
      }else{
        return(TRUE) #pair passed check
      }
    }		
  }
  
  all_pairs <- pair_matrix(1:dim(gene)[2])
  all_pair_names <- pair_matrix(colnames(gene))
  
  if(verbose){cat("There are", nrow(all_pairs), "possible marker pairs to test.\n")}
  
  check_one_pair <- function(pair_num){
    pair <- all_pairs[pair_num,]
    pass_checks <- check_linkage(m1 = gene[,pair[1]], m2 = gene[,pair[2]], thresh_param = thresh_param)
    return(pass_checks)
  }	
  
  check_multi_pairs <- function(pair_V){
    pair_checks <- unlist(lapply(pair_V, check_one_pair))
    return(pair_checks)
  }
  
  #chunk up the jobs based on how many cores we want to use
  pair_list <- chunkV(1:nrow(all_pairs), n_cores)
  
  if (run_parallel) {
    
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    cape_dir_full <- find.package("cape")
    cape_dir <- gsub("cape_pkg/cape","cape_pkg", cape_dir_full)
    clusterExport(cl, "cape_dir", envir=environment())
    clusterEvalQ(cl, .libPaths(cape_dir))
    good_pair_list <- foreach(p = 1:length(pair_list)) %dopar% {
      check_multi_pairs(pair_list[[p]])
    }
    stopCluster(cl)
    
  } else {
    
    good_pair_list <- c()
    index <- 1:length(pair_list)
    for (p in index) {
      good_pair_list <- rbind(good_pair_list, check_multi_pairs(pair_list[[p]]))
    }
  }
  
  testV <- unlist(good_pair_list)
  idxV <- unlist(pair_list)
  
  pairs_mat <- all_pair_names[idxV[which(testV)],,drop = FALSE]
  if(verbose){
    cat(dim(pairs_mat)[1], "marker pairs will be tested.\n")
  }
  
  colnames(pairs_mat) <- c("marker1", "marker2")
  rownames(pairs_mat) <- NULL
  pairs_mat <- rbind(pairs_mat, all_covar_pairs)
  
  if(verbose){
    cat(dim(pairs_mat)[1], "pairs including covariates will be tested.\n")
  }
  
  return(pairs_mat)
  
}