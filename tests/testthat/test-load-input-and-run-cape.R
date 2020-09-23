context("test load input and run cape")
skip_on_cran()

# Code for running cape with a yaml string
# We first read the parameter file found in demo/demo_qtl and then
# we read the file into a string which is then passed to the 
# load_input_and_run_cape function.

#===============================================================
# load all the necessary libraries
#===============================================================
needed.packages <- c("evd", "Matrix", "fdrtool", "shape", "corpcor", "RColorBrewer", "doParallel",
                     "foreach", "caTools", "stringr", "abind", "propagate", "here", "testthat")

for(i in 1:length(needed.packages)){library(needed.packages[i], character.only = TRUE)}

# load the yaml file as a string
demo_path <- here("demo/demo_qtl")
param_file <- file.path(demo_path, "NON_NZO.parameters.yml")
yaml_parameters <- readLines(param_file)

results_path <- here("tests/testthat/results")

# before running cape, make sure the results directory is clean
unlink(list.files(results_path))
data_path <- here("tests/testthat/testdata/demo_qtl_data")

# we get the datafile
data_file <- file.path(data_path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")

final_cross <- load_input_and_run_cape(input_file = data_file, yaml_params = yaml_parameters, 
                                       results_path = results_path)

test_that("test that the output is not null", {
  expect_true(!is.null(final_cross))  # not null
})
