set.seed(1234)

library(cape)

results_path <- here::here("demo", "demo_qtl")
data_path <- here::here("tests", "testthat", "testdata", "demo_qtl_data")

data_file <- file.path(data_path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")
param_file <- file.path(results_path, "NON_NZO.parameters.yml")

cross <- read_population(data_file)
cross_obj <- cape2mpp(cross)
data_obj <- cross_obj$data_obj
geno_obj <- cross_obj$geno_obj$geno

final_cross <- run_cape(data_obj, geno_obj, p_or_q = 0.05, verbose = TRUE, 
                param_file = param_file, results_path = results_path)

plot_variant_influences(final_cross)
