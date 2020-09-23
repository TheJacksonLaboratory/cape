context("test get pairs for pairscan")
library(here)
library(doParallel)

# test that requires sourcing the files under test
pkg_path <- here("R")
for (f in list.files(pkg_path)) {
  source(file.path(pkg_path, f))
}

#This example finds marker pairs to test in a randomly
#generated genotype matrix
possible_genotypes <- c(0, 0.5, 1)
genotype_matrix <- matrix(sample(possible_genotypes, 1000, replace = TRUE), nrow = 100, ncol = 10)
colnames(genotype_matrix) <- paste("marker", letters[1:10], sep = "_")
pairs_which <- get_pairs_for_pairscan(genotype_matrix, min_per_genotype = 6, verbose = TRUE)

test_that("test that there are no self pairs", {
  len <- length(pairs_which[pairs_which[, "marker1"] == pairs_which[, "marker2"]])
  expect_equal(0, len)
})

test_that("test that the result set size is correct", {
  expect_lt(10, dim(pairs_which)[1])
  expect_equal(2, dim(pairs_which)[2])
})