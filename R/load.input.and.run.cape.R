#'
#' This function loads the input file path and runs cape
#'
#' @param input_file data input to be loaded 
#' @param yaml_params a parameter set up in the form of a YAML string
#' @param results_path path to the results
#' @param use_kinship TRUE or FALSE, FALSE by default
#' @param run_parallel boolean, if TRUE runs certain parts of the code as parallel blocks
#' @param results_file the name of the saved data.obj RData file. The base name is used as the base name for all saved RData files.
#' @param p_or_q A threshold indicating the maximum adjusted p value considered 
#' @param n_cores integer, default is 4
#' @param run_singlescan boolean, defaul: TRUE
#' @param run_pairscan boolean, default: TRUE
#' @param error_prop_coef, boolean, default: TRUE
#' @param error_prop_perm, boolean, default: TRUE
#' @param initialize_only, boolean, default: FALSE
#' @param verbose boolean, output goes to stdout
#'
#' @export
load.input.and.run.cape <- function(input_file = NULL, yaml_params = NULL, results_path = NULL,
                                    use_kinship = FALSE, run_parallel = FALSE, results_file = "cross.RData", p_or_q = 0.05, 
                                    n_cores = 4, run_singlescan = TRUE, run_pairscan = TRUE, error_prop_coef = TRUE, 
                                    error_prop_perm = TRUE, initialize_only = FALSE, verbose = TRUE){
  		
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
  
  final.cross <- run.cape(data.obj, geno.obj, results.file = results_file, p.or.q = p_or_q, 
                          n.cores = n_cores, run.singlescan = run_singlescan, run.pairscan = run_pairscan, 
                          error.prop.coef = error_prop_coef, error.prop.perm = error_prop_perm, 
                          initialize.only = initialize_only, verbose = verbose, run.parallel = run_parallel,
                          yaml.params = yaml_params, results.paths = results_path)
}
