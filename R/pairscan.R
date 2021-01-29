#' This function performs the pairwise scan on all markers.
#'
#' This function performs the pairwise regression on all selected marker pairs.
#' The phenotypes used can be either eigentraits or raw phenotypes. Permutation
#' testing is also performed.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' @param scan_what A character string uniquely identifying whether eigentraits
#'   or raw traits should be scanned. Options are "eigentraits", "raw_traits"
#' @param pairscan_null_size The total size of the null distribution.
#' This is DIFFERENT than the number of permutations to run. Each permutation
#' generates n choose 2 elements for the pairscan. So for example, a permutation
#' that tests 100 pairs of markers will generate a null distribution of size 4950.
#' This process is repeated until the total null size is reached. If the null size
#' is set to 5000, two permutations of 100 markers would be done to get to a null
#' distribution size of 5000.
#' @param max_pair_cor A numeric value between 0 and 1 indicating the maximum
#'   Pearson correlation that two markers are allowed. If the correlation
#'   between a pair of markers exceeds this threshold, the pair is not tested.
#'   If this value is set to NULL, min_per_genotype must have a numeric value.
#' @param min_per_genotype The minimum number of individuals allowable per
#'   genotype combination. If for a given marker pair, one of the genotype combinations is
#'   underrepresented, the marker pair is not tested. If this value is NULL,
#'   max_pair_cor must have a numeric value.
#' @param kin_obj a kinship object calculated by \code{\link{kinship}}.
#' @param num_pairs_limit A number indicating the maximum number of pairs to
#'   scan. If the number of pairs exceeds this threshold, the function asks for
#'   confirmation before proceeding with the pairwise scan.
#' @param num_perm_limit A number indicating the maximum number of total
#'   permutations that will be performed. If the number of total permutations
#'   exceeds this threshold, the function asks for confirmation before
#'   proceeding with the pairwise scan.
#' @param overwrite_alert If TRUE raises a warning to users not to overwrite 
#'   their data object with a singlescan object. A warning necessary after a 
#'   new version of cape began separating results from different functions into
#'   different results objects
#' @param run_parallel Whether to run the analysis on parallel CPUs
#' @param n_cores The number of CPUs to use if run_parallel is TRUE
#' @param verbose Whether to write progress to the screen
#'
#'
#' @details Not all marker pairs are necessarily tested. Before markers are
#'   tested for interaction, they are checked for several conditions. Pairs are
#'   discarded if (1) at least one of the markers is on the X chromosome, or (2)
#'   there are fewer than min_per_genotype individuals in any of the genotype
#'   combinations.
#'
#' @return This function returns an object assigned to pairscan_obj in 
#' \code{\link{run_cape}}.
#'
#' The results object is a list of five elements:
#' ref_allele: The allele used as the reference for the tests.
#' max_pair_cor: The maximum pairwise correlation between marker pairs
#' pairscan_results: A list with one element per trait. The element for
#' each trait is a list of the following three elements:
#'    pairscan_effects: the effect sizes from the linear models
#'    pairscan_se: the standard errors from the linear models
#'    model_covariance: the model covariance from the linear models.
#' pairscan_perm: The same structure as pairscan_results, but for the
#' permuted data.
#' pairs_tested_perm: A matrix of the marker pairs used in the permutation
#' tests.
#'   
#' @seealso \code{\link{select_markers_for_pairscan}}, \code{\link{plot_pairscan}}
#' 
#' @examples 
#' \dontrun{
#' pairscan_obj <- pairscan(data_obj, geno_obj, pairscan_null_size = 10000)
#' }
#'
#' @export
pairscan <- function(data_obj, geno_obj = NULL,
  scan_what = c("eigentraits", "raw_traits"), pairscan_null_size = NULL, 
  max_pair_cor = NULL, min_per_genotype = NULL, kin_obj = NULL, 
  num_pairs_limit = 1e6, num_perm_limit = 1e7, overwrite_alert = TRUE, 
  run_parallel = FALSE, n_cores = 4, verbose = FALSE) {
  
  marker_selection_method <- data_obj$marker_selection_method
  
  if(!run_parallel){n_cores = 1}
  
  use_kinship <- data_obj$use_kinship
  
  
  if(overwrite_alert){
    choice <- readline(prompt = "Please make sure you assign the output 
                       of this function to a pairscan_obj, and NOT the data_obj. It will 
                       overwrite the data_obj.\nDo you want to continue (y/n) ")
    if(choice == "n"){stop()}
  }
  
  
  pairscan_obj <- list()
  pairscan_obj$ref_allele <- data_obj$ref_allele
  pairscan_obj$max_pair_cor <- data_obj$max_pair_cor
  pairscan_obj$min_per_genotype <- data_obj$min_per_genotype
  
  if(is.null(pairscan_null_size)){
    stop("The final size of the null distribution must be specified.")
  }
  
  pheno <- get_pheno(data_obj, scan_what)	
  
  covar_info <- get_covar(data_obj)
  covar_names <- covar_info$covar_names
  covar_table <- covar_info$covar_table
  
  #find the phenotypic covariates. These will
  #be tested separately, and not as part of a
  #chromosome
  
  if(is.null(data_obj$geno_for_pairscan)){
    stop("select_markers_for_pairscan() must be run before pairscan()")
  }
  
  
  #add the covariates (geno and pheno) 
  #to the genotype matrix so that we 
  #test all pairs
  gene <- get_geno_with_covar(data_obj, geno_obj, g_covar = TRUE, p_covar = TRUE, 
    for_pairscan = TRUE)	
  
  #fill in a matrix to index the marker pairs
  if(verbose){cat("Getting marker pairs for pairscan...\n")}
  pared_marker_mat <- get_pairs_for_pairscan(
    gene,
    covar_names,
    max_pair_cor,
    min_per_genotype,
    run_parallel = run_parallel,
    n_cores = n_cores,
    verbose = verbose
  )
  
  num_pairs <- dim(pared_marker_mat)[1]
  
  if(num_pairs == 0){
    stop("There are no pairs to test. Try raising max_pair_cor or reducing 
         min_per_genotype.")
  }
  
  if(!is.null(num_pairs_limit) && num_pairs > num_pairs_limit){
    warning("\nThe number of pairs (",num_pairs,") exceeds ", num_pairs_limit, ".\n", sep = "")
    go_on <- readline(prompt = "Do you want to continue (y/n)?\n")
    if(length(grep("n", go_on))){
      message("Stopping pairwise scan...\n")
      return(pairscan_obj)
    }else{
      message("Continuing pairwise scan...\n")
    }
  }
  
  if(verbose){cat("Performing pairwise tests...\n")}
  #run one_pairscan for each phenotype with results in scanone_result
  if(!use_kinship){
    pairscan_results <- pairscan_noKin(data_obj, pheno_mat = pheno, 
      geno_mat = gene, covar_table = covar_table, marker_pairs = pared_marker_mat, 
      n_perm = pairscan_null_size, verbose = verbose, n_cores = n_cores, 
      run_parallel = run_parallel)
  }else{
    pairscan_results <- pairscan_kin(data_obj, geno_obj = geno_obj, 
      scan_what = scan_what, marker_pairs = pared_marker_mat, kin_obj = kin_obj, 
      verbose = verbose, run_parallel = run_parallel, n_cores = n_cores)
  }	
  
  # print(str(pairscan_results))
  
  pairscan_obj$pairscan_results <- pairscan_results	
  
  if(pairscan_null_size > 0){	
    if(use_kinship){
      pairscan_perm <- pairscan_null_kin(data_obj, geno_obj, kin_obj, 
        scan_what = scan_what, pairscan_null_size = pairscan_null_size, 
        max_pair_cor = max_pair_cor, min_per_genotype, verbose = verbose, 
        marker_selection_method = marker_selection_method, 
        run_parallel = run_parallel, n_cores = n_cores)
    }else{
      pairscan_perm <- pairscan_null(data_obj, geno_obj, scan_what = scan_what, 
        pairscan_null_size = pairscan_null_size, max_pair_cor = max_pair_cor, 
        min_per_genotype, verbose = verbose, marker_selection_method = marker_selection_method, 
        run_parallel = run_parallel, n_cores = n_cores)
    }
    #add the results to the data object
    pairscan_obj$pairscan_perm <- pairscan_perm$pairscan_perm 
    #add the results to the data object
    pairscan_obj$pairs_tested_perm <- pairscan_perm$pairs_tested_perm 
    
  }
  
  return(pairscan_obj) #and return it
  
  }
