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
#' @param param_file path to yml parameter file for running cape
#' @param create_report boolean, if true we create the corresponding HTML report page
#' @param qtl_id_col argument for read_population, an optional column number for individual IDs
#' @param qtl_na_strings argument for read_population, an optional string for missing values
#' 
#' @import here
#' @importFrom qtl2 read_cross2
#' 
#' @export
load_input_and_run_cape <- function(input_file = NULL, yaml_params = NULL, results_path = NULL,
                                    run_parallel = FALSE, results_file = "cross.RDS", p_or_q = 0.05, 
                                    n_cores = 4, initialize_only = FALSE, verbose = TRUE, param_file = NULL,
                                    create_report = FALSE, qtl_id_col = NULL, qtl_na_strings = "-"){
  # create log file
  cape_log <- file(file.path(results_path, "cape.log"), open="wt")
  sink(cape_log)
  sink(cape_log, type="message")
  
  if (endsWith(input_file, ".yaml") || endsWith(input_file, ".json") || endsWith(input_file, ".yml")) {
    # QTL2 file type (with json/yml control file in a folder)
    cat("Read_cross2: Read file into a qtl2 object\n")
    qtl2 <- read_cross2(input_file)
    cat("Transform QTL2 object to Cape object\n")
    cape_object <- qtl2_to_cape(qtl2)
    data_obj <- cape_object$data_obj
    geno_obj <- cape_object$geno_obj 
  } else if (endsWith(input_file, ".csv")){
    # QTL file type as a single CSV file
    cat("Read population into cross object\n")
    cross <- read_population(input_file, id_col = qtl_id_col, 
	    na_strings = qtl_na_strings, verbose = verbose)
    cat("Convert from Cape1 object to Cape2 object\n")
    cross_obj <- cape2mpp(cross)
    data_obj <- cross_obj$data_obj
    geno_obj <- cross_obj$geno_obj$geno
  } else if (endsWith(input_file, ".RDS") || endsWith(input_file, ".rds") 
             || endsWith(input_file, ".RDATA") || endsWith(input_file, ".rdata")) {
    datafile <- list.files(input_file, pattern="data", full.names=TRUE)
    if(datafile == "") {
      datafile <- list.files(input_file, pattern="pheno")
    }
    cat(paste("Data file is: ", datafile, "\n"))

    genofile <- list.files(input_file, pattern="geno", full.names=TRUE)
    if (genofile == "") {
      stop("No Geno file (file with a \'geno\' in its name) was found.")
    }
    cat(paste("Geno file is: ", genofile, "\n"))
    
    cat("Reading data file...\n")
    cape_obj <- readRDS(datafile)
    cat("Reading Geno file...\n")
    geno_obj <- readRDS(genofile)
    
    # genotype coding
    het_val <- 0.3 #could be 0.5, but I think we can go as low as 0.3
    dom_geno <- geno_obj
    for(i in 1:dim(dom_geno)[3]){
      dom_mat <- dom_geno[,,i]
      dom_mat[which(dom_mat >= het_val)] <- 1
      dom_mat[which(dom_mat < het_val)] <- 0
      dom_geno[,,i] <- dom_mat
    }
    cat("Creating CAPE object\n")
    data_obj <- Cape$new(
      parameter_file = param_file,
      yaml_parameters = yaml_params,
      results_path = results_path,
      pheno = cape_obj$pheno,
      chromosome = cape_obj$chromosome,
      marker_num = cape_obj$marker_num,
      marker_location = cape_obj$marker_location,
      geno_names = cape_obj$geno_names,
      geno = dom_geno
    )
    
  } else if (endsWith(input_file, ".ped")) {
    # file is PLINK data type and we load the corresponding ped, map and pheno files
    # ped = input_file

    # cross_obj <- plink2cape(ped, map, pheno, out, overwrite = TRUE)

    # data_obj <- cross_obj$data_obj
    # geno_obj <- cross_obj$geno_obj$geno
    stop("PLINK data is not yet supported")
  }
  
  final_cross <- run_cape(data_obj, geno_obj, results_file = results_file, p_or_q = p_or_q, 
                          n_cores = n_cores, initialize_only = initialize_only, verbose = verbose, 
                          run_parallel = run_parallel, param_file = param_file,
                          yaml_params = yaml_params, results_path = results_path, plot_pdf = FALSE)
  
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
