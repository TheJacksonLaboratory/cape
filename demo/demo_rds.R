set.seed(1234)

library(cape)

results_path <- here::here("demo", "demo_rds")
data_path <- here::here("tests", "testthat", "testdata", "demo_rds_data")

data_file <- file.path(data_path, "cape_data.RDS")
geno_file <- file.path(data_path, "cape_geno.RDS")

cape_obj <- readRDS(data_file)
cape_geno <- readRDS(geno_file)

param_file <- file.path(results_path, "NON_NZO.parameters.yml")

# genotype coding
het_val <- 0.3 #could be 0.5, but I think we can go as low as 0.3
dom_geno <- cape_geno
for(i in 1:dim(dom_geno)[3]){
  dom_mat <- dom_geno[,,i]
  dom_mat[which(dom_mat >= het_val)] <- 1
  dom_mat[which(dom_mat < het_val)] <- 0
  dom_geno[,,i] <- dom_mat
}

cross_obj <- Cape$new(
  parameter_file = param_file,
  results_path = results_path,
  pheno = cape_obj$pheno,
  chromosome = cape_obj$chromosome,
  marker_num = cape_obj$marker_num,
  marker_location = cape_obj$marker_location,
  geno_names = cape_obj$geno_names,
  geno = dom_geno
)

final_cross <- run_cape(pheno_obj = cross_obj, geno_obj = cape_geno, p_or_q = 0.05,
                        results_path = results_path, param_file = param_file)

plot_variant_influences(final_cross)
