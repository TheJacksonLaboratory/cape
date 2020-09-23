context("test CAPE class")
library(here)

pkg_path <- here("R")
for (f in list.files(pkg_path)) {
  source(file.path(pkg_path, f))
}

test_data_path <- here("tests/testthat/testdata/demo_qtl_data")
params_data_path <- here("demo/demo_qtl")

file_name <- file.path(test_data_path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")
param_file <- file.path(params_data_path, "NON_NZO.parameters.yml")
results_path <- file.path("results")
dir.create(results_path, showWarnings = FALSE)

cross <- read_population(file_name)
cross_obj <- cape2mpp(cross)
geno_obj <- cross_obj$geno_obj$geno

data_obj <- Cape$new(
  parameter_file = param_file,
  results_path = results_path,
  pheno = cross_obj$data_obj$pheno,
  chromosome = cross_obj$data_obj$chromosome,
  marker_num = cross_obj$data_obj$marker_num,
  marker_location = cross_obj$data_obj$marker_location,
  geno_names = dimnames(geno_obj),
  geno = geno_obj
)

test_that("test that the pheno matrix is right", {
  expect_equal(208, dim(data_obj$pheno)[1])
})

# change the reference allele, and check the change
test_that("check setting the ref_allele", {
  data_obj$ref_allele <- "B"
  expect_equal("B", data_obj$ref_allele)
})

# check that setting the wrong type throws and error
# test_that("throw error on bad type", {
#   expect_error(cp$ref_allele <- c(1,2,3))
# })


