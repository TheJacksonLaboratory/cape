context("test CAPE class")

test.data.path <- here("tests/testthat/testdata")
file.name <- file.path(test.data.path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")
param.file <- file.path(test.data.path, "cape.parameters.yml")
results.path <- file.path("results")
dir.create(results.path, showWarnings = FALSE)

cross <- read.population(file.name)
cross.obj <- cape2mpp(cross)
geno.obj <- cross.obj$geno.obj$geno

data.obj <- Cape$new(
  parameter_file = param.file,
  results_path = results.path,
  pheno = cross.obj$data.obj$pheno,
  chromosome = cross.obj$data.obj$chromosome,
  marker_num = cross.obj$data.obj$marker_num,
  marker_location = cross.obj$data.obj$marker_location,
  geno_names = dimnames(geno.obj),
  geno = geno.obj
)

test_that("test that the pheno matrix is right", {
  expect_equal(208, dim(data.obj$pheno)[1])
})

# change the reference allele, and check the change
test_that("check setting the ref_allele", {
  data.obj$ref_allele <- "B"
  expect_equal("B", data.obj$ref_allele)
})

# check that setting the wrong type throws and error
# test_that("throw error on bad type", {
#   expect_error(cp$ref_allele <- c(1,2,3))
# })


