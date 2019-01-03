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

library(yaml)
library(R6)
library(microbenchmark)

test.data.path <- "/Users/emersj/projects/cape/cape/tests/testthat/testdata"
file.name <- file.path(test.data.path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")
param.file <- file.path(test.data.path, "cape.parameters.yml")

cross <- read.population(file.name)
cross.obj <- cape2mpp(cross)

cape.obj <- new(
  "Cape",
  pheno = cross.obj$data.obj$pheno,
  chromosome = cross.obj$data.obj$chromosome,
  marker_num = cross.obj$data.obj$marker.num,
  marker_location = cross.obj$data.obj$marker.location,
  geno_names = cross.obj$data.obj$geno.names,
  ref_allele = "A",
  geno = cross.obj$geno.obj$geno,
  parameters = yaml.load_file(param.file),
  covar_table = array(),
  flat_geno = array(),
  non_allelic_covar = c(NA_character_)
  
)

tracemem(cape.obj)

c6 <- CapeR6$new(
  pheno = cross.obj$data.obj$pheno,
  geno = cross.obj$geno.obj$geno,
  covar_table = array()
)

print(
microbenchmark(
  cpCreateCovarTable(cape.obj) <- c("LEP_24", "total_fat"),
  c6$cteate_covar_table(c("LEP_24", "total_fat")),
  times=5
)
)




