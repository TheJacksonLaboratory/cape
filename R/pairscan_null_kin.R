#' Generates a null distribution for the pairscan
#' 
#' This function generates a null distribution
#' for the pairscan. For each permutation,
#' it runs a single scan and selects markers
#' in the same manner as for the true test.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' @param kin_obj a kinship object
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
#' @param min_per_geno The minimum number of individuals allowable per
#'   genotype. If for a given marker pair, one of the genotypes is
#'   underrepresented, the marker pair is not tested. If this value is NULL,
#'   max_pair_cor must have a numeric value.
#' @param model_family Indicates the model family of the phenotypes. This can be 
#'   either "gaussian" or "binomial".
#' @param marker_selection_method options are "top_effects", "uniform", "effects_dist", "by_gene"
#' @param run_parallel Whether to run the analysis on multiple CPUs.
#' @param n_cores The number of CPUs to use if run_parallel is TRUE.
#' @param verbose Whether to print progress to the screen. Defaults to FALSE.
#' 
#' @return This function returns a list with two elements, one containing
#' the results of the permutations, and the other containing the markers
#' that were tested in the individual permutations.
#' @keywords internal
#' 
pairscan_null_kin <- function(data_obj, geno_obj = NULL, kin_obj = NULL, 
  scan_what = c("eigentraits", "raw_traits"), pairscan_null_size = NULL, 
  max_pair_cor = NULL, min_per_geno = NULL, model_family = "gaussian", 
  marker_selection_method = c("top_effects", "uniform", "effects_dist", "by_gene"), 
  run_parallel = FALSE, n_cores = 4, verbose = FALSE){
  
  marker_selection_method <- data_obj$marker_selection_method
  
  ref_allele <- data_obj$ref_allele
  
  
  if(is.null(pairscan_null_size)){
    stop("The total number of permutations must be specified.")
  }
  
  
  #If the user does not specify a scan_what, 
  #default to eigentraits, basically, if eigen,
  #et, or ET are anywhere in the string, use the
  #eigentraits, otherwise, use raw phenotypes
  type_choice <- c(grep("eig", scan_what, ignore.case = TRUE), grep("ET", scan_what, ignore.case = TRUE)) #look for any version of eigen or eigentrait, the user might use.
  if(length(type_choice) > 0){ #if we find any, use the eigentrait matrix
    pheno <- data_obj$ET
  }else{
    pheno <- data_obj$pheno #otherwise, use the raw phenotype matrix
  }
  
  num_pheno <- dim(pheno)[2]
  #use the full genotype matrix to select 
  #markers for generating the null in the 
  #pairscan
  
  geno <- get_geno(data_obj, geno_obj)
  
  
  #make a list to hold the results. 
  #Tables from each of the phenotypes will be
  #put into this list
  results_perm_list <- vector(mode = "list", length = num_pheno)
  
  #generate the null distribution for the pairscan 
  #do a singlescan on each permuted trait.
  #find the top n markers
  #combine these into a unique list
  #do the pairscan on these markers for all traits
  
  if(verbose){
    cat("\nGenerating null distribution...\n")
  }
  
  all_pairs_tested <- NULL
  
  n_top_markers <- ncol(data_obj$geno_for_pairscan)
  final_perm <- 1
  while(final_perm < pairscan_null_size){
    
    # TODO what is happening here? should this be a deep copy?
    perm_data_obj <- data_obj
    perm_data_obj$pheno <- perm_data_obj$pheno[sample(nrow(perm_data_obj$pheno)),]
    
    if(marker_selection_method != "by_gene" && marker_selection_method != "from_list"){
      single_scan_result <- array(NA, dim = c(length(data_obj$geno_names[[3]]), num_pheno, (dim(geno)[[2]]-1)))
      dimnames(single_scan_result) <- list(data_obj$geno_names[[3]], colnames(pheno), dimnames(geno)[[2]][-which(dimnames(geno)[[2]] == ref_allele)])
      
      if(verbose){cat("Performing single marker scans of permuted traits.\n")}
      
      #sink all the warnings about solutions close to zero to a file
      #one_singlescan <- singlescan(perm_data_obj, geno_obj, kin_obj, n_perm = 0, 
      #  model_family = model_family, run_parallel = run_parallel, n_cores = n_cores, 
      #  verbose = verbose, overwrite_alert = FALSE)

      one_singlescan <- singlescan(perm_data_obj, geno_obj, kin_obj = NULL, n_perm = 0, 
        model_family = model_family, run_parallel = run_parallel, n_cores = n_cores, 
        verbose = verbose, overwrite_alert = FALSE)
      

      single_scan_result <- one_singlescan$singlescan_t_stats
      
      if(verbose){cat("Selecting markers for permuted pairscan...\n")}				
      #use this singlescan to select markers for a permuted pairscan
      
      if(marker_selection_method == "top_effects"){
        perm_data_obj <- select_markers_for_pairscan(perm_data_obj, singlescan_obj = single_scan_result, geno_obj, 
          num_alleles = n_top_markers, peak_density = data_obj$peak_density, window_size = data_obj$window_size, 
          tolerance = data_obj$tolerance, plot_peaks = FALSE, verbose = verbose)
      }
    }else{ 
      
      # if(marker_selection_method == "by_gene"){
        # #if we are using a gene-based method
        # #use a permuted gene list to select
        # #SNPs near genes
        # perm_data_obj <- select_markers_for_pairscan.by_gene(perm_data_obj, ref_allele = ref_allele, geno_obj = geno_obj, 
                                                             # gene.list = sample(gene.list), num.snps = ncol(data_obj$geno_for_pairscan), 
                                                             # organism = data_obj$organism)
      # }
      if(marker_selection_method == "from_list"){
        single_scan_result <- list("ref_allele" = ref_allele)
        specific_markers <- colnames(perm_data_obj$geno_for_pairscan)
        perm_data_obj <- select_markers_for_pairscan(data_obj, singlescan_obj = single_scan_result, geno_obj, specific_markers = specific_markers)
      }
    }
    
    
    if(verbose){cat("\tGetting markers for permuted pairscan...\n")}
    top_marker_pairs <- get_pairs_for_pairscan(gene = perm_data_obj$geno_for_pairscan, 
    max_pair_cor = max_pair_cor, min_per_genotype = min_per_geno, 
    run_parallel = run_parallel, n_cores = n_cores, verbose = verbose)
    total_pairs <- nrow(top_marker_pairs)
    num_to_add <- 10
    #we don't want to do more permutations than specified
    #so trim the final pair matrix down to get only
    #the specified number of permutations plus a few
    #because some pairs are always rejected
    
    if(final_perm+dim(top_marker_pairs)[1] > pairscan_null_size){
      num_needed <- pairscan_null_size - final_perm
      #testing just one pair was messing this up, so 
      #always test at least two pairs
      top_marker_pairs <- top_marker_pairs[1:(num_needed+(min(c(num_to_add, total_pairs)))),,drop=FALSE]
    }
    
    if(verbose){cat("\tTesting", dim(top_marker_pairs)[1], "pairs...\n")}
    all_pairs_tested <- rbind(all_pairs_tested, top_marker_pairs)
    
    #run the pairscan for each permuted phenotype and the pairs we just found
    if(verbose){cat("Performing marker pair scans of permuted traits with kinship correction...\n")}
    
    
    #run a pairscan on these markers and each permuted phenotype
    pairscan_results <- pairscan_kin(perm_data_obj, geno_obj, scan_what = scan_what, 
    marker_pairs = top_marker_pairs, kin_obj, run_parallel = run_parallel, 
    n_cores = n_cores, verbose = verbose)
    
    #integrate the results into the permutation object
    one_perm <- pairscan_results[[1]]
    #because there will be different numbers of markers each time, just take 
    #the marker names, the intercept, and the effects for marker1 marker2 and 
    #their interaction
    # last.col = dim(one_perm[[1]])[2]
    # take.col <- c(1:3, (last.col-2):last.col)
    if(final_perm == 1){ #if this is the first time through, 
      #just copy the results into the results_perm_list
      for(p in 1:num_pheno){
        results_perm_list[[p]] <- pairscan_results[[p]]
      }
    }else{
      if(!is.null(one_perm)){
        for(p in 1:num_pheno){
          for(i in 1:length(pairscan_results[[1]])){
            results_perm_list[[p]][[i]] <- rbind(results_perm_list[[p]][[i]], pairscan_results[[p]][[i]])
          }
        }
      }
    }
    
    final_perm <- nrow(results_perm_list[[1]][[1]]) #end looping through phenotypes
    if(verbose){cat("\t", final_perm, " null tests: ", round((final_perm/pairscan_null_size)*100), "%...\n", sep = "")} 
  } #end when we have enough permutations
  
  
  names(results_perm_list) <- colnames(pheno)
  results_list <- list("pairscan_perm" = results_perm_list, "pairs_tested_perm" = all_pairs_tested)
  return(results_list)
  
}
