#'
#' This function instantiates an R6 Cape object and is meant to be called from Python
#' using RPy2
#'
#' @param yaml.parameters a parameter set up in the form of a YAML string
#' @param yaml.parameters path to the results
#' @param pheno phenotype data
#' @param chromosome chromosomes list
#' @param marker.num marker data
#' @param marker.location marker location data
#' @param geno.names genotype data object
#' @param geno genotype data
#' @param use.kinship TRUE or FALSE, FALSE by default
#'
#' @return \codeCape object
#'
#' @export
instantiate.cape <- function(yaml.parameters = NULL, results_path = NULL, pheno = NULL, chromosome = NULL, 
                             marker.num = NULL, marker.location = NULL, geno.names = NULL, geno = NULL, use.kinship = FALSE){
  data.obj <- Cape$new(
    yaml_parameters = yaml.parameters,
    results_path = here("results"),
    pheno = pheno,
    chromosome = chromosome,
    marker_num = marker.num,
    marker_location = marker.location,
    geno_names = geno.names,
    geno = geno,
    use_kinship = use.kinship
  )
  return $data.obj
  
}