#' Loads input and run CAPE
#' 
#' This function loads the input file path and runs cape
#' It is used to run CAPE from a non R script (python)
#'
#' @param input_file data input to be loaded 
#' @param yaml_params a parameter set up in the form of a YAML string
#' @param results_path path to the results
#' @param run_parallel boolean, if TRUE runs certain parts of the code as parallel blocks
#' @param results_file the name of the saved data_obj RDS file. The base name is used as the base name for all saved RDS files.
#' @param p_or_q A threshold indicating the maximum adjusted p value considered 
#' @param n_cores integer, default is 4
#' @param initialize_only boolean, default: FALSE
#' @param verbose boolean, output goes to stdout
#' @param create_report boolean, if true we create the corresponding HTML report page
#' 
#' @import here
#' @importFrom qtl2 read_cross2
#' 
#' @examples 
#' \dontrun{
#' #load input in qtl2 zip format
#' load_input_and_run_cape("cross_file.zip")
#' 
#' #load input in qtl csv format
#' load_input_and_run_cape("cross_file.csv")
#' }
#' 
#' @export
load_input_and_run_cape <- function(input_file = NULL, yaml_params = NULL, results_path = NULL,
                                    run_parallel = FALSE, results_file = "cross.RDS", p_or_q = 0.05, 
                                    n_cores = 4, initialize_only = FALSE, verbose = TRUE, create_report = FALSE){

  if (endsWith(input_file, ".yaml") || endsWith(input_file, ".json") || endsWith(input_file, ".yml")) {
    # QTL2 file type (with json/yml control file in a folder)
    qtl2 <- read_cross2(input_file)

    cape_object <- qtl2_to_cape(qtl2)
    data_obj <- cape_object$data_obj
    geno_obj <- cape_object$geno_obj 
  } else if (endsWith(input_file, ".csv")){
    # QTL file type as a single CSV file
    cross <- read_population(input_file)
    cross_obj <- cape2mpp(cross)
    data_obj <- cross_obj$data_obj
    geno_obj <- cross_obj$geno_obj$geno
  } else if (endsWith(input_file, ".ped")) {
    # file is PLINK data type and we load the corresponding ped, map and pheno files
    # ped = input_file

    # cross_obj <- plink2cape(ped, map, pheno, out, overwrite = TRUE)

    # data_obj <- cross_obj$data_obj
    # geno_obj <- cross_obj$geno_obj$geno
    stop("PLINK data is not yet supported")
  }

  # create log file
  cape_log <- file(file.path(results_path, "cape.log"), open="wt")
  sink(cape_log)
  sink(cape_log, type="message")
  
  final_cross <- run_cape(data_obj, geno_obj, results_file = results_file, p_or_q = p_or_q, 
                          n_cores = n_cores, initialize_only = initialize_only, verbose = verbose, run_parallel = run_parallel,
                          yaml_params = yaml_params, results_path = results_path)
  
  if(create_report) {
    cat("Rendering result page...\n")
    # copy result page rmd to result folder. The file in in the resource folder of the project running this method.
    # we hardcoded the place where the markdown file is located on the production VM
    file.copy("/opt/cape/cape/cape_results.Rmd", results_path, overwrite = TRUE)
    # render result page
    rmarkdown::render(file.path(results_path, "cape_results.Rmd"), params = list(results_dir = results_path))
    cat("Result page rendered.\n")
  }
  
  cat("The CAPE analysis has been successfully completed.")
  ## reset message sink and close the file connection
  sink(type="message")
  sink()
  close(cape_log)
  
  return(final_cross)
}
