# # where are the data?
# plink.path <- "/Users/emersj/projects/cape/data/plink"
# 
# # which file?
# file.base <- "test*"
# 
# # find the `map` and `ped` files in the globbed directory6
# map <- NULL
# ped <- NULL
# 
# for (file.name in Sys.glob(file.path(plink.path, file.base))) {
#   file.name <- tolower(file.name)
#   if ("ped" == file_ext(file.name)) {
#     ped <- file.name
#   } else if ("map" == file_ext(file.name)) {
#     bed <- file.name
#   }
# }
# 
# # drop out if the files are missing
# if (is.null(bed) | is.null(ped)) {
#   stop("Input path is missing the BED or PED file.")
# }
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




