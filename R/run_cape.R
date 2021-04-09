#' Runs CAPE
#'
#' This function takes in a data object and genotype object that
#' have been formatted for cape, as well as a string identifying
#' a parameter file. It runs cape on the data using the parameters
#' specified in the file.
#' 
#' This function assumes you already have all required libraries 
#' and functions loaded.
#'
#' @param pheno_obj the cape object holding the phenotype data returned by 
#' @param geno_obj the genotype object
#' @param results_file the name of the saved data_obj RDS file. The base 
#' name is used as the base name for all saved RDS files.
#' @param p_or_q A threshold indicating the maximum adjusted p value 
#' considered significant. Or, if FDR p value correction is used, the
#' the maximum q value considered significant.
#' @param n_cores The number of CPUs to use if run_parallel is set to TRUE
#' @param initialize_only, If TRUE, cape will not be run. Instead an initialized
#' data object will be returned. This data object will contain normalized and mean-centered
#' trait values, and eigentraits, and will have covariates specified. However, the 
#' singlescan, pairscan, etc. will not be run.
#' @param verbose Whether progress should be printed to the screen
#' @param run_parallel boolean, if TRUE runs certain parts of the code as parallel blocks
#' @param param_file yaml full path to the parameter file
#' @param yaml_params yaml string containing the parameters. Either the param_file or 
#' yaml_params can be null.
#' @param results_path path that results should be written to.
#'
#' @return This function invisibly returns the data object with all final 
#' data included. In addition, data saved to the data_obj$results_path directory
#' 
#' @importFrom utils read.table
#' 
#' @examples 
#' \dontrun{
#' final_data_obj <- run_cape(pheno_obj, geno_obj)
#' }
#' 
#' @export
run_cape <- function(pheno_obj, geno_obj, 
  results_file = "cross.RDS", p_or_q = 0.05, n_cores = 4, 
  initialize_only = FALSE, verbose = TRUE, run_parallel = FALSE, 
  param_file = NULL, yaml_params = NULL, results_path = NULL){
  
  # Instantiate the Cape R6 object
  data_obj <- Cape$new(
  		parameter_file = param_file,
  		yaml_parameters = yaml_params,
  		results_path = results_path,
  		pheno = pheno_obj$pheno,
  		chromosome = pheno_obj$chromosome,
  		marker_num = pheno_obj$marker_num,
  		marker_location = pheno_obj$marker_location,
  		geno_names = pheno_obj$geno_names,
  		geno = geno_obj
  )
  
  results_base_name <- gsub(".RDS", "", results_file)
    
  # since this is the main data_obj, we can't allow it to return FALSE, 
  #check for the file first
  prior_data_obj <- data_obj$read_rds(results_file)
  if (isFALSE(prior_data_obj)) {
    data_obj <- compare_markers(data_obj, geno_obj)
  } else {
    # things can get pretty confusing if these values don't match between 
    #the parameter file and the old data_obj
    prior_data_obj$save_results <- data_obj$save_results
    prior_data_obj$use_saved_results <- data_obj$use_saved_results
    data_obj <- prior_data_obj
  }
  
  #===============================================================
  # figure out how to synchronize get_eigentraits and
  # the calculation of the kinship matrix. They both may
  # remove individuals. Maybe trim the kinship object to match
  # the eigentraits
  #===============================================================
  #if we want to use a kinship correction
  if(as.logical(data_obj$use_kinship)){
    kin_file_name <- paste0(results_base_name, "_kinship.RDS")
    kin_obj <- data_obj$read_rds(kin_file_name)
    
    if (isFALSE(kin_obj)){
      #if there isn't a kinship object already, we need to make one
      kin_obj <- kinship(data_obj, geno_obj, type = "overall", 
      pop = data_obj$pop)
      data_obj$save_rds(kin_obj, kin_file_name)
    }

    #===============================================================
    # We need a complete genotype matrix to calculate the kinship
    # adjusted genotypes later on.
    # Check for missing values in the genotype matrix.
    # If there are missing values, impute them.
    # Write out the imputed matrix, or read this in if it already
    # exists.
    #===============================================================
    #we need to impute the missing values
    imp_data_file <- paste0(results_base_name, "_data_imputed.RDS")
    imp_geno_file <- paste0(results_base_name, "_geno_imputed.RDS")

    # check if there is already a saved genotype object
    geno <- data_obj$read_rds(imp_geno_file)

    if (isFALSE(geno)) {  #if the imputation hasn't been done already
      geno <- get_geno(data_obj, geno_obj)
      missing_vals <- which(is.na(geno))

      if (length(missing_vals) > 0) { #if there are missing values, impute them
        message("There are missing values in geno_obj. Running impute_missing_geno...\n")
        geno_imp <- impute_missing_geno(data_obj, geno_obj = geno_obj, k = 10, 
        	ind_missing_thresh = 0, marker_missing_thresh = 0, prioritize = "fewer",
        	max_region_size = NULL, min_region_size = NULL, run_parallel = run_parallel,
        	verbose = verbose, n_cores = n_cores)

        # update and save the data_obj
        data_obj <- geno_imp$data_obj
        data_obj$save_rds(data_obj, imp_data_file)

        # update and save the geno_obj
        geno_obj <- geno_imp$geno_obj
        data_obj$save_rds(geno_obj, imp_geno_file)
      
        # recalculate the kinship matrix with the updated objects
        #kin_obj <- kinship(data_obj, geno_obj, type = "overall", pop = data_obj$pop)
        #data_obj$save_rds(kin_obj, kin_file_name)

      } #end case for when there are missing values but no imputed genotypes

    } else { #if the imputation has been done, then it must have been done for the data_obj too
      #data_obj <- geno_imp$data_obj
      geno_obj <- geno
    }

  }
    
  if(verbose){cat("Removing unused markers...\n")}
  data_obj <- remove_unused_markers(data_obj, geno_obj, verbose = verbose)
  combined_data_obj <- delete_underscore(data_obj, geno_obj)

  data_obj <- combined_data_obj$data_obj
  geno_obj <- combined_data_obj$geno_obj
  
  #because the genotype object can be changed by the above step, 
  #save the final version. (or change the above step so it doesn't 
  final_geno_file <- paste0(results_base_name, "_geno.RDS")
  data_obj$save_rds(geno_obj, final_geno_file)

  #str(data_obj$geno_names)
  #str(dimnames(geno_obj))

  if(!is.null(data_obj$covariates)){
    data_obj <- pheno2covar(data_obj, data_obj$covariates)
  }
  if(!is.null(data_obj$marker_covariates)){
    data_obj <- marker2covar(data_obj, geno_obj, markers = data_obj$marker_covariates)
  }
  
  data_obj <- select_pheno(data_obj, pheno_which = data_obj$traits)	
  
  if(length(grep("eig", data_obj$scan_what, ignore.case = TRUE)) > 0){
    data_obj <- get_eigentraits(
      data_obj, 
      scale_pheno = as.logical(data_obj$traits_scaled), 
      normalize_pheno = as.logical(data_obj$traits_normalized)
    )
    
    if(data_obj$save_results){
      data_obj$plotSVD("svd.pdf")
      data_obj$plotSVD("svd.jpg")
    }
    
    # TODO update select_eigentraits
    data_obj <- select_eigentraits(data_obj, traits_which = data_obj$eig_which)
  }
  
  data_obj$save_rds(data_obj, results_file)

  if(initialize_only){
    return(data_obj)
  }
  
  #===============================================================
  # run singlescan
  #===============================================================
  singlescan_results_file <- paste0(results_base_name, "_singlescan.RDS")
  
  singlescan_obj <- data_obj$read_rds(singlescan_results_file)
  
  if(!data_obj$use_kinship) {
    kin_obj <- NULL
  }
  
  if (isFALSE(singlescan_obj)) {
      
      singlescan_obj <- singlescan(
        data_obj, geno_obj, kin_obj = kin_obj, n_perm = data_obj$singlescan_perm,
        alpha = data_obj$alpha, verbose = verbose, run_parallel = run_parallel,
        n_cores = n_cores, model_family = "gaussian", overwrite_alert = FALSE
      )
      
      data_obj$save_rds(singlescan_obj, singlescan_results_file)
      
      if(data_obj$save_results){
        for(ph in 1:ncol(singlescan_obj$singlescan_effects)){
          filename <- paste0("Singlescan_", colnames(singlescan_obj$singlescan_effects)[ph], "_Standardized.jpg")
          data_obj$plotSinglescan(filename, singlescan_obj, width = 20, height = 6, 
            units = "in", res = 300, standardized = TRUE, allele_labels = NULL, 
            alpha = data_obj$alpha, include_covars = TRUE, line_type = "l", pch = 16, cex = 0.5, 
            lwd = 3, traits = colnames(singlescan_obj$singlescan_effects)[ph])
        }
      }
      
      if(data_obj$save_results){
        for(ph in 1:ncol(singlescan_obj$singlescan_effects)){
          filename <- paste0("Singlescan_", colnames(singlescan_obj$singlescan_effects)[ph], "_Effects.jpg")
          data_obj$plotSinglescan(filename, singlescan_obj, width = 20, height = 6, units = "in", res = 300,
            standardized = FALSE, allele_labels = NULL, alpha = data_obj$alpha, include_covars = TRUE,
            line_type = "l", pch = 16, cex = 0.5, lwd = 3, traits = colnames(singlescan_obj$singlescan_effects)[ph])
        }
      }
  }
  
  #===============================================================
  # run pairscan
  #===============================================================
  pairscan_file <- paste0(results_base_name, "_pairscan.RDS")
  
  pairscan_obj <- data_obj$read_rds(pairscan_file)
  
  if (isFALSE(pairscan_obj) | is.null(data_obj$geno_for_pairscan)) {
      
      marker_selection_method <- data_obj$marker_selection_method
      num_alleles_in_pairscan <- data_obj$num_alleles_in_pairscan
      peak_density <- data_obj$peak_density
      max_pair_cor <- data_obj$max_pair_cor
      min_per_genotype <- data_obj$min_per_genotype
      pairscan_null_size <- data_obj$pairscan_null_size
      scan_what <- data_obj$scan_what
      if(marker_selection_method == "top_effects"){
        data_obj <- select_markers_for_pairscan(data_obj, singlescan_obj, geno_obj,
          num_alleles = num_alleles_in_pairscan, peak_density = peak_density,
          verbose = verbose, plot_peaks = FALSE)
      }
      
      if(marker_selection_method == "from_list"){
        if(is.null(data_obj$snp_file)){stop("Please specify a marker list in the parameter file.\n")}
        snp_file <- file.path(results_path, data_obj$snp_file)
        if(!file.exists(snp_file)){stop("Can't fine the specified marker list.\n")}
        specific_markers <- as.matrix(read.table(snp_file, sep = "\t", stringsAsFactors = FALSE))
        data_obj <- select_markers_for_pairscan(data_obj, singlescan_obj, geno_obj,
          specific_markers = specific_markers[,1], verbose = verbose, plot_peaks = FALSE)
      }
      
      # if(marker_selection_method == "by_gene"){
        # gene_list_mat <- read.table("gene_list.txt", sep = "\t", stringsAsFactors = FALSE)		
        # gene_list <- gene_list_mat[,1]
        # data_obj <- select_markers_for_pairscan_by_gene(data_obj, geno_obj, gene_list = gene_list, 
                                                        # bp_buffer = data_obj$bp_buffer, organism = data_obj$organism)
      # } else {
        # gene_list <- NULL
      # }
      
      data_obj$save_rds(data_obj, results_file)
      
      pairscan_obj <- pairscan(data_obj, geno_obj, scan_what = scan_what, 
        pairscan_null_size = pairscan_null_size, min_per_genotype = min_per_genotype, 
        max_pair_cor = max_pair_cor, verbose = verbose, num_pairs_limit = Inf, 
        overwrite_alert = FALSE, run_parallel = run_parallel, n_cores = n_cores, 
        kin_obj = kin_obj)
      
      data_obj$save_rds(pairscan_obj, pairscan_file)
      
      if(data_obj$save_results){
        data_obj$plotPairscan("Pairscan_Regression.pdf", pairscan_obj, 
          phenotype = NULL, show_marker_labels = TRUE, show_alleles = FALSE)
        data_obj$plotPairscan("Pairscan_Regression.jpg", pairscan_obj, 
          phenotype = NULL, show_marker_labels = TRUE, show_alleles = FALSE)
      }

      data_obj$save_rds(data_obj, results_file)

  }
  
  #===============================================================
  # run reparametrization
  #===============================================================
  

  if(!data_obj$use_saved_results || is.null(data_obj$var_to_var_influences)){
    data_obj <- error_prop(data_obj, pairscan_obj, perm = FALSE, verbose = verbose,
      n_cores = n_cores, run_parallel = run_parallel)
    data_obj$save_rds(data_obj, results_file)
  }
  
  if(!data_obj$use_saved_results || is.null(data_obj$var_to_var_influences_perm)){	
    data_obj <- error_prop(data_obj, pairscan_obj, perm = TRUE, verbose = verbose,
      n_cores = n_cores, run_parallel = run_parallel)
     data_obj$save_rds(data_obj, results_file)
  }
  
  if(!data_obj$use_saved_results || is.null(data_obj$var_to_var_p_val)){
    data_obj <- calc_p(data_obj, pval_correction = data_obj$pval_correction)
  }
  
  #if(length(grep("e", data_obj$scan_what, ignore_case = TRUE)) > 0){
  #  transform_to_phenospace <- TRUE
  #}else{
  #  transform_to_phenospace <- FALSE	
  #}
  
  if(!data_obj$use_saved_results || is.null(data_obj$max_var_to_pheno_influence)){
    data_obj <- direct_influence(data_obj, pairscan_obj, 
      transform_to_phenospace = data_obj$transform_to_phenospace, verbose = verbose, 
      pval_correction = data_obj$pval_correction, save_permutations = TRUE, 
      n_cores = n_cores)
      data_obj$save_rds(data_obj, results_file)
  }
  
  if(data_obj$save_results){
    data_obj$writeVariantInfluences("Variant_Influences.csv", p_or_q = max(c(p_or_q, 0.2)))

    data_obj$writeVariantInfluences("Variant_Influences_Interactions.csv", 
      include_main_effects = FALSE, p_or_q = max(c(p_or_q, 0.2)))
  
    data_obj$plotVariantInfluences("variant_influences.pdf", width = 10, height = 7,
      p_or_q = p_or_q, standardize = FALSE, not_tested_col = "lightgray", 
      covar_width = NULL, pheno_width = NULL)

    data_obj$plotVariantInfluences("variant_influences.jpg", width = 10, height = 7,
      p_or_q = p_or_q, standardize = FALSE, not_tested_col = "lightgray", 
      covar_width = NULL, pheno_width = NULL)
  }

  if(!data_obj$use_saved_results || is.null(data_obj$full_net)){
    data_obj <- get_network(data_obj, geno_obj, p_or_q = p_or_q, 
    collapse_linked_markers = FALSE)
  }
  
  if(!data_obj$use_saved_results || is.null(data_obj$collapsed_net)){
    data_obj <- get_network(data_obj, geno_obj, p_or_q = p_or_q, threshold_power = 1, 
    collapse_linked_markers = TRUE, plot_linkage_blocks = FALSE)
  }
  
  data_obj$save_rds(data_obj, results_file)
  
  if(data_obj$save_results){
    data_obj$plotNetwork("Network_Circular.pdf", label_gap = 10, label_cex = 1.5, show_alleles = FALSE)
    data_obj$plotNetwork("Network_Circular.jpg", label_gap = 10, label_cex = 1.5, show_alleles = FALSE)
  
    if(dim(geno_obj)[2] == 8){
      data_obj$plotNetwork("Network_Circular_DO.pdf", label_gap = 10, label_cex = 1.5, show_alleles = TRUE)
      data_obj$plotNetwork("Network_Circular_DO.jpg", label_gap = 10, label_cex = 1.5, show_alleles = TRUE)
    }	
  
    data_obj$plotFullNetwork("Network_View.pdf", zoom = 1.2, node_radius = 0.3, 
      label_nodes = TRUE, label_offset = 0.4, label_cex = 0.5, bg_col = "lightgray", 
      arrow_length = 0.1, layout_matrix = "layout_with_kk", legend_position = "topright", 
      edge_lwd = 1, legend_radius = 2, legend_cex = 0.7, xshift = -1)
  
    data_obj$plotFullNetwork("Network_View.jpg", zoom = 1.2, node_radius = 0.3, 
      label_nodes = TRUE, label_offset = 0.4, label_cex = 0.5, bg_col = "lightgray", 
      arrow_length = 0.1, layout_matrix = "layout_with_kk", legend_position = "topright", 
      edge_lwd = 1, legend_radius = 2, legend_cex = 0.7, xshift = -1)
  }
  
  data_obj$save_rds(data_obj, results_file)
  
  invisible(data_obj)
  
  
}
