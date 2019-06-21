# # where are the data?
# plink.path <- "/Users/emersj/projects/cape/data/plink"
# 
# # which file?
# file.base <- "test*"
# 
# # find the `map` and `ped` files in the globbed directory
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


test.data.path <- "/Users/emersj/projects/cape/cape/tests/testthat/testdata"
file.name <- file.path(test.data.path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")
param.file <- file.path(test.data.path, "cape.parameters.yml")

cross <- read.population(file.name)
cross.obj <- cape2mpp(cross)
geno.obj <- cross.obj$geno.obj$geno

data.obj <- Cape$new(
  parameter_file = param.file,
  results_path = "/Users/emersj/projects/cape/cape/results",
  pheno = cross.obj$data.obj$pheno,
  chromosome = cross.obj$data.obj$chromosome,
  marker_num = cross.obj$data.obj$marker_num,
  marker_location = cross.obj$data.obj$marker_location,
  geno_names = dimnames(geno.obj),
  geno = geno.obj
)

# TODO search throughout the code for data.obj[[number]]
# TODO that call won't work because the data.obj is an environmnet not a list now

# TODO search for all the variables in read.parameters and prepend with "data.obj$"
# TODO because these are part of the cape object now, and not in the global namespace

# TODO handle all the calls to .RData files with readRDS and come up 
# TODO with a scheme for controlling their storage and retrival via a cape$method()

# TODO remove all calls to require() and ensure that the libraries are in the DESCRIPTION file

final.cross <- run.cape(data.obj, geno.obj, p.or.q = 0.05, path = ".", results.file = "cross.RData",
                        n.cores = 4, run.singlescan = TRUE, run.pairscan = TRUE, error.prop.coef = TRUE,
                        error.prop.perm = TRUE, initialize.only = FALSE, verbose = TRUE, run.parallel = TRUE)




