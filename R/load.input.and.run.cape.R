#'
#' This function loads the input file path and runs cape
#'
#' @param input_file data input to be loaded 
#' @param yaml_params a parameter set up in the form of a YAML string
#' @param results_path path to the results
#' @param use_kinship TRUE or FALSE, FALSE by default
#' @param run_parallel FALSE by default
#'
#'
#' @export
load.input.and.run.cape <- function(input_file = NULL, yaml_params = NULL, results_path = NULL, use_kinship = FALSE,
                                    run.parallel = FALSE){
  if(!require(here)){install.packages("here")}
  
  test.data.path <- here("tests/testthat/testdata")
  # "NON_NZO_Reifsnyder_pgm_CAPE_num.csv"
  # file.name <- file.path(test.data.path, input_file)
  file.name <- file.path(input_file) 
  
  # param.file <- file.path(test.data.path, "cape.parameters.yml")
  
  cross <- read.population(file.name)
  cross.obj <- cape2mpp(cross)
  geno.obj <- cross.obj$geno.obj$geno
  
  data.obj <- Cape$new(
    # parameter_file = param.file,
    yaml_parameters = yaml_params,
    results_path = results_path,
    pheno = cross.obj$data.obj$pheno,
    chromosome = cross.obj$data.obj$chromosome,
    marker_num = cross.obj$data.obj$marker_num,
    marker_location = cross.obj$data.obj$marker_location,
    geno_names = dimnames(geno.obj),
    geno = geno.obj,
    use_kinship = use_kinship
  )
  
  # TODO remove all calls to require() and ensure that the libraries are in the DESCRIPTION file
  
  snp.file = here("tests/testthat/testdata/NON_NZO_marker_list.txt")
  
  final.cross <- run.cape(data.obj, geno.obj, results.file = "cross.RData", p.or.q = 0.05, snp.file = NULL,
                          n.cores = 4, run.singlescan = TRUE, run.pairscan = TRUE, error.prop.coef = TRUE,
                          error.prop.perm = TRUE, initialize.only = FALSE, verbose = TRUE, run.parallel = run_parallel)
  
}
