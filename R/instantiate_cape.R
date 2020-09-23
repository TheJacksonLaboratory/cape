#'
#' This function instantiates an R6 Cape object and is meant to be called from Python
#' using RPy2
#'
#' @param param_file YAML parameter file 
#' @param yaml_params a parameter set up in the form of a YAML string
#' @param results_path path to the results
#' @param pheno phenotype data
#' @param chromosome chromosomes list
#' @param marker_num marker data
#' @param marker_location marker location data
#' @param geno_names genotype data object
#' @param geno genotype data
#' @param use_kinship TRUE or FALSE, FALSE by default
#'
#' @return Cape object
#'
#' @export
instantiate_cape <- function(param_file = NULL, yaml_params = NULL, results_path = NULL, pheno = NULL, chromosome = NULL, 
                             marker_num = NULL, marker_location = NULL, geno_names = NULL, geno = NULL, use_kinship = FALSE){
  data_obj <- Cape$new(
    parameter_file = param_file,
    yaml_parameters = yaml_params,
    results_path = here("results"),
    pheno = pheno,
    chromosome = chromosome,
    marker_num = marker_num,
    marker_location = marker_location,
    geno_names = geno_names,
    geno = geno,
    use_kinship = use_kinship
  )
  return $data_obj
  
}