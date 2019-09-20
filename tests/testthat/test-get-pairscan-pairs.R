context("test get pairs for pairscan")
library(here)
library(doParallel)

# test that requires sourcing the files under test
pkg.path <- here("R")
for (f in list.files(pkg.path)) {
  source(file.path(pkg.path, f))
}

#This example finds marker pairs to test in a randomly
#generated genotype matrix
possible.genotypes <- c(0, 0.5, 1)
genotype.matrix <- matrix(sample(possible.genotypes, 1000, replace = TRUE), nrow = 100, ncol = 10)
colnames(genotype.matrix) <- paste("marker", letters[1:10], sep = "_")
pairs.which <- get.pairs.for.pairscan(genotype.matrix, min.per.genotype = 6, verbose = TRUE)

test_that("test that there are no self pairs", {
  len <- length(pairs.which[pairs.which[, "marker1"] == pairs.which[, "marker2"]])
  expect_equal(0, len)
})

# change the reference allele, and check the change
test_that("test that the result set size is correct", {
  expect_lt(20, dim(pairs.which)[1])
  expect_equal(2, dim(pairs.which)[2])
})