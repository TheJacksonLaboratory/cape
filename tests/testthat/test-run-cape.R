# context("test run cape")
# skip_on_cran()
#
# #code for running cape
# #make sure there is a file called cape.parameters.txt in the results directory
#
# #===============================================================
# # load all the necessary libraries
# #===============================================================
# needed.packages <- c("evd", "Matrix", "fdrtool", "shape", "corpcor", "RColorBrewer", "doParallel",
#                      "foreach", "caTools", "stringr", "abind", "propagate", "here", "testthat")
#
# for(i in 1:length(needed.packages)){library(needed.packages[i], character.only = TRUE)}
#
#
# data_dir <- here("tests", "testthat", "testdata")
#
# #===============================================================
# # set up directories
# #===============================================================
# testing = FALSE #if this is set to TRUE, downsample to 1000 markers for testing
# cape_parameter_file <- file.path(data_dir, "cape.parameters.txt")
#
# cape_do_dir <- here("capempp")
# results_dir <- file.path(data_dir, "ResultsTest")
# n_cores <- detectCores() - 1
# run_parallel <- TRUE
# #===============================================================
#
#
# #===============================================================
# #source cape code
# #===============================================================
# # all_fun <- list.files(pattern = ".R", path = cape_do_dir, full.names = TRUE)
# # for(i in 1:length(all.fun)){source(all_fun[i])}
#
# #===============================================================
# # read in the data
# #===============================================================
#
# cross <- read_population(file.path(data_dir, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv"))
# cross_obj <- cape2mpp(cross)
# cross <- cross_obj$data_obj
# geno <- cross_obj$geno_obj
#
#
# #===============================================================
# #set the results directory and run cape from the parameter file
# #===============================================================
# # before running cape, make sure the results directory is clean
# unlink(list.files(results_dir))
#
# final_cross <- run_cape(cross, geno, param_file = cape_parameter_file, p_or_q = 0.05, n_cores = n_cores,
#                         error_prop_coef = TRUE, error_prop_perm = TRUE, verbose = TRUE)
#
#
# test_that("test that the output is not null", {
#   expect_true(!is.null(final_cross))  # not null
# })
