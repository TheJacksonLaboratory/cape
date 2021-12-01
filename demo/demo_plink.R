set.seed(1234)

library(cape)

results_path <- here::here("demo", "demo_PLINK")
data_path <- here::here("tests", "testthat", "testdata", "demo_PLINK_data")

ped <- file.path(data_path, "test.ped")
map <- file.path(data_path, "test.map")
pheno <- file.path(data_path, "test.pheno")
out <- file.path(data_path, "test.csv")
param_file <- file.path(results_path, "plink.parameters.yml")

cross_obj <- plink2cape(ped, map, pheno, out, overwrite = TRUE)

data_obj <- cross_obj$data_obj
geno_obj <- cross_obj$geno_obj$geno

final_cross <- run_cape(pheno_obj = data_obj, geno_obj, 
	p_or_q = 0.05, verbose = TRUE, param_file = param_file, 
	results_path = results_path)

plot_variant_influences(final_cross)
