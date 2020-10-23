#' Loads input and run CAPE
#' 
#' This function loads the input file path and runs cape
#' It is used to run CAPE from a non R script (python)
#'
#' @param input_file data input to be loaded 
#' @param yaml_params a parameter set up in the form of a YAML string
#' @param results_path path to the results
#' @param run_parallel boolean, if TRUE runs certain parts of the code as parallel blocks
#' @param results_file the name of the saved data_obj RData file. The base name is used as the base name for all saved RData files.
#' @param p_or_q A threshold indicating the maximum adjusted p value considered 
#' @param n_cores integer, default is 4
#' @param initialize_only boolean, default: FALSE
#' @param verbose boolean, output goes to stdout
#' 
#' @import here
#' 
#' @export
load_input_and_run_cape <- function(input_file = NULL, yaml_params = NULL, results_path = NULL,
                                    run_parallel = FALSE, results_file = "cross.RData", p_or_q = 0.05, 
                                    n_cores = 4, initialize_only = FALSE, verbose = TRUE){
  		
  # if(!require(here)){install.packages("here")}
  
  # if R/QTL2 file format
  if (endsWith(input_file, ".zip")) {

    qtl2 <- read_cross2(input_file)
    
    cape_object <- qtl2_to_cape(qtl2)
    data_obj <- cape_object$data_obj
    geno_obj <- cape_object$geno_obj 
  } else if (endsWith(input_file, ".csv")){
    # csv file format like NON_NZO...csv
    cross <- read_population(input_file)
    cross_obj <- cape2mpp(cross)
	  data_obj <- cross_obj$data_obj
    geno_obj <- cross_obj$geno_obj$geno
  }

  final_cross <- run_cape(data_obj, geno_obj, results_file = results_file, p_or_q = p_or_q, 
                          n_cores = n_cores, initialize_only = initialize_only, verbose = verbose, run_parallel = run_parallel,
                          yaml_params = yaml_params, results_path = results_path)
}
