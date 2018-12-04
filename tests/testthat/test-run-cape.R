context("test run cape")

#code for running cape
#make sure there is a file called cape.parameters.txt in the results directory

#===============================================================
# load all the necessary libraries
#===============================================================
needed.packages <- c("evd", "Matrix", "fdrtool", "shape", "corpcor", "RColorBrewer", "doParallel",
                     "foreach", "caTools", "stringr", "abind", "propagate", "here", "testthat")

for(i in 1:length(needed.packages)){library(needed.packages[i], character.only = TRUE)}


data.dir <- here("tests", "testthat", "testdata")

#===============================================================
# set up directories
#===============================================================
testing = FALSE #if this is set to TRUE, downsample to 1000 markers for testing
cape.parameter.file <- file.path(data.dir, "cape.parameters.txt")

cape.do.dir <- here("capempp")
results.dir <- file.path(data.dir, "ResultsTest")
n.cores <- detectCores() - 1
run.parallel <- TRUE
#===============================================================


#===============================================================
#source cape code
#===============================================================
all.fun <- list.files(pattern = ".R", path = cape.do.dir, full.names = TRUE)
for(i in 1:length(all.fun)){source(all.fun[i])}

#===============================================================
# read in the data
#===============================================================

cross <- read.population(file.path(data.dir, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv"))
cross.obj <- cape2mpp(cross)
cross <- cross.obj$data.obj
geno <- cross.obj$geno.obj


#===============================================================
#set the results directory and run cape from the parameter file
#===============================================================
# before running cape, make sure the results directory is clean
unlink(list.files(results.dir))

setwd(results.dir)

final.cross <- run.cape(cross, geno, cape.parameter.file, p.or.q = 0.05, n.cores = n.cores, 
                        run.singlescan = TRUE, run.pairscan = TRUE, error.prop.coef = TRUE, 
                        error.prop.perm = TRUE, verbose = TRUE)


test_that("test that the output is not null", {
  expect_true(!is.null(final.cross))  # not null
})
