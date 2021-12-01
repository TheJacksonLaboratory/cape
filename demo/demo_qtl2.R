set.seed(1234)

library(cape)

demo.path <- here::here("demo", "demo_qtl2")
data_path <- here::here("tests", "testthat", "testdata", "demo_qtl2_data")

data_file <- file.path(data_path, "iron.yaml")
param_file <- file.path(demo.path, "iron.parameters.yml")

#read in example qtl2 data from remote host
iron.qtl2 <- read_cross2(data_file)

iron.cape <- qtl2_to_cape(cross = iron.qtl2)
data_obj <- iron.cape$data_obj
geno_obj <- iron.cape$geno_obj 

final_cross <- run_cape(pheno_obj = data_obj, geno_obj, 
	p_or_q = 0.05, verbose = TRUE, param_file = param_file, 
	results_path = demo.path)

plot_variant_influences(final_cross)
