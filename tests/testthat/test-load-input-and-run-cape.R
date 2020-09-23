context("test load input and run cape")
skip_on_cran()

#code for running cape with a yaml string
#make sure there is a file called cape.parameters.txt in the results directory

#===============================================================
# load all the necessary libraries
#===============================================================
needed.packages <- c("evd", "Matrix", "fdrtool", "shape", "corpcor", "RColorBrewer", "doParallel",
                     "foreach", "caTools", "stringr", "abind", "propagate", "here", "testthat")

for(i in 1:length(needed.packages)){library(needed.packages[i], character.only = TRUE)}


yaml_parameters <- "#================================================
# General Parameters 
#================================================
traits:
  - BW_24
  - INS_24
  - log_GLU_24
covariates:
  - pgm
scan_what:
  - eigentraits 
traits_scaled:
  - true
traits_normalized:
  - true
eig_which:
  - 1
- 2 
pval_correction:
  - fdr 
use_kinship:
  - false
pop:
  - 2PP
save_results:
  - true
use_saved_results:
  - false
transform_to_phenospace:
  - true

#================================================
# Single Scan Parameters 
#================================================
ref_allele:
  - A 
singlescan_perm:
  - 0
alpha:
  - 0.05
- 0.01

#================================================
# Marker Selection Parameters 
#================================================
marker_selection_method:
  - top_effects
peak_density:
  - 0.5 
tolerance:
  - 5 
num_alleles_in_pairscan:
  - 84
#================================================
# Pairscan Parameters 
#================================================
max_pair_cor:
  - 0.5 
pairscan_null_size:
  - 5000
"

#results_path <- file.path("results")
results_path <- here("tests/testthat/results")

# before running cape, make sure the results directory is clean
unlink(list.files(results_path))
data_path <- here("tests/testthat/testdata/demo_qtl_data")

data_file <- file.path(data_path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")

final_cross <- load_input_and_run_cape(input_file = data_file, yaml_params = yaml_parameters, results_path = results_path,
                                       run_parallel = FALSE, results_file = "cross.RData", p_or_q = 0.05, 
                                       n_cores = 4, initialize_only = FALSE, verbose = TRUE)


test_that("test that the output is not null", {
  expect_true(!is.null(final_cross))  # not null
})
