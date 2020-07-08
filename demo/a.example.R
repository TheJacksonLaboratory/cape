if(!require(here)){install.packages("here")}
# if(!require(cape)){install.packages("cape")}

test.data.path <- here("tests/testthat/testdata")
file.name <- file.path(test.data.path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")
param.file <- file.path(test.data.path, "cape.parameters.yml")

cross <- read.population(file.name)
cross.obj <- cape2mpp(cross)
data.obj <- cross.obj$data.obj
geno.obj <- cross.obj$geno.obj$geno

# TODO remove all calls to require() and ensure that the libraries are in the DESCRIPTION file

snp.file = here("tests/testthat/testdata/NON_NZO_marker_list.txt")

final.cross <- run.cape(data.obj, geno.obj, results.file = "cross.RData", p.or.q = 0.05, snp.file = NULL,
                        n.cores = 4, run.singlescan = TRUE, run.pairscan = TRUE, error.prop.coef = TRUE,
                        error.prop.perm = TRUE, initialize.only = FALSE, verbose = TRUE, run.parallel = FALSE,
                        param.file = param.file, results.path = here("results"))
