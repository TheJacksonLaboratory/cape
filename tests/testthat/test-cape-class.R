context("test CAPE class")
library(here)

# instantiate the class
arr <- array(data     = 1:27,
             dim      = c(3, 3, 3),
             dimnames = list(c("a", "b", "c"),
                             c("d", "e", "f"),
                             c("g", "h", "i")))

cp <- new(
  "Cape",
  pheno = matrix(1:9, nrow = 3, ncol = 3),
  chromosome = c("1", "1", "1"),
  marker_num = as.integer(c(1, 2, 3)),
  marker_location = c(1.1, 2.2, 3.3),
  geno_names = list("A", "B"),
  ref_allele = "A",
  geno = arr,
  parameters = arr
)

test_that("test that the pheno matrix is right", {
  expect_equal(3, dim(getPheno(cp))[1])
})

# change the reference allele, and check the change
test_that("check setting the ref_allele", {
  setRefAllele(cp) <- "B"
  expect_equal("B", getRefAllele(cp))
})

# check that setting the wrong type throws and error
test_that("throw error on bad type", {
  expect_error(setRefAllele(cp) <- c(1,2,3))
})


