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
if(!require(here)){install.packages("here")}

library(cape, lib.loc = "/opt/cape/cape_pkg")

test.data.path <- here("tests/testthat/testdata")
# file.name <- file.path(test.data.path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")
param.file <- file.path(test.data.path, "cape.parameters.chesler.DO.exp.group.yml")
# param.file <- "/Users/emersj/projects/cape/cape/tests/testthat/testdata/cape.parameters.chesler.DO.exp.group.yml"

# cross <- read.population(file.name)
# cross.obj <- cape2mpp(cross)
# geno.obj <- cross.obj$geno.obj$geno

pheno.obj<-readRDS(here("../data/CheslerDO_Data/CheslerDO.pheno.RData")) 
geno.obj<-readRDS(here("../data/CheslerDO_Data/CheslerDO.geno.RData"))

# set all the dotted list names in the pheno.obj to snake case
names(pheno.obj) <- gsub("[.]", "_", names(pheno.obj))

# combined.data.obj <- cape::delete.underscore(pheno.obj, geno.obj)
# 
# pheno.obj <- combined.data.obj$data.obj
# geno.obj <- combined.data.obj$geno.obj

#grep("X", data.obj$chromosome, ignore.case = TRUE)

pheno.obj <- cape::remove.unused.markers(pheno.obj, geno.obj)

browser()

data.obj <- Cape$new(
  parameter_file = param.file,
  results_path = here("results"),
  pheno = pheno.obj$pheno,
  chromosome = pheno.obj$chromosome,
  marker_num = pheno.obj$marker_num,
  marker_location = pheno.obj$marker_location,
  geno_names = pheno.obj$geno_names,
  geno = geno.obj,
  use_kinship = TRUE
)

# TODO remove all calls to require() and ensure that the libraries are in the DESCRIPTION file


final.cross <- run.cape(data.obj, geno.obj, results.file = "cross.RData", p.or.q = 0.05, snp.file = NULL,
                        n.cores = 4, run.singlescan = TRUE, run.pairscan = TRUE, error.prop.coef = TRUE,
                        error.prop.perm = TRUE, initialize.only = FALSE, verbose = TRUE, run.parallel = TRUE)




