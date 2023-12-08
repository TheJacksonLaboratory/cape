#' Read the parameter file, add missing entries
#' 
#' This function returns reads in the YAML file
#' and checks for any parameters that might not
#' be included. This may not matter for the given
#' run, but it's handy to be able to check for 
#' any and all potential variables.
#'
#' @param filename full path to the .yml file holding CAPE parameters (is not needed if yaml_parameters is provided)
#' @param yaml_parameters yaml string holding CAPE parameters (can be NULL)
#'
#' @return Returns a named list with all possible options
#' 
#' @import yaml
#' 
#' @export
read_parameters <- function(filename = "cape.parameters.yml", yaml_parameters = NULL){
  if (!is.null(yaml_parameters)) {
    parameter_table <- yaml.load(yaml_parameters)
  } else {
    parameter_table <- read_yaml(filename)
  }
  
  #================================================
  # general parameters
  #================================================
  gen_param <- c("traits", "covariates", "marker_covariates", "traits_scaled", "traits_normalized", "scan_what", "eig_which", "use_kinship", "kinship_type", "locus", "pop", "pval_correction", "transform_to_phenospace")
  
  #================================================
  #single scan parameters
  #================================================
  single_param <- c("ref_allele", "singlescan_perm", "alpha")
  
  #================================================
  # marker selection
  #================================================
  marker.param <- c("marker_selection_method", "snp_file", "peak_density", "tolerance", "window_size", "num_alleles_in_pairscan", "bp_buffer", "organism")
  
  #================================================
  # pair scan
  #================================================
  pair_param <- c("max_pair_cor", "min_per_genotype", "pairscan_null_size")
  
  for ( param_name in c(gen_param, single_param, marker.param, pair_param) ) {
    if ( !is.element(param_name, names(parameter_table)) ){
      parameter_table[param_name] = NULL
    }
  }
  
  return(parameter_table)
  
}
