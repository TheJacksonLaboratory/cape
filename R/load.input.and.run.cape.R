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
load.input.and.run.cape <- function(input_file = NULL, yaml_params = NULL, results_path = NULL,
	use_kinship = FALSE, run_parallel = FALSE, results.file = "cross.RData", p.or.q = 0.05, 
  	n.cores = 4, run.singlescan = TRUE, run.pairscan = TRUE, error.prop.coef = TRUE, 
  	error.prop.perm = TRUE, initialize.only = FALSE, verbose = TRUE){
  		
  if(!require(here)){install.packages("here")}
  
  # if R/QTL2 file format
  if (endsWith(input_file, ".zip")) {

    qtl2 <- read_cross2(input_file)
    
    cape.object <- qtl2_to_cape(qtl2)
    data.obj <- cape.object$data.obj
    geno.obj <- cape.object$geno.obj 
  } else if (endsWith(input_file, ".csv")){
    # csv file format like NON_NZO...csv
    cross <- read.population(input_file)
    cross.obj <- cape2mpp(cross)
	data.obj <- cross.obj$data.obj
    geno.obj <- cross.obj$geno.obj$geno
  }
  
  # TODO remove all calls to require() and ensure that the libraries are in the DESCRIPTION file
  final.cross <- run.cape(data.obj, geno.obj, results.file = "cross.RData", p.or.q = 0.05, 
  	n.cores = 4, run.singlescan = run.singlescan, run.pairscan = run.pairscan, 
  	error.prop.coef = error.prop.coef, error.prop.perm = error.prop.perm, 
  	initialize.only = initialize.only, verbose = verbose, run.parallel = run_parallel,
  	yaml.params = yaml_params, results.paths = results_path)
  
}
