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

pheno <- cross$pheno
chromosome <- cross$chromosome
marker_num <- cross$marker.num
marker_location <- cross$marker.location
geno_names <- cross$marker.names
ref_allele <- "A"  # TODO this should come from the parameter file
geno <- cross$geno
parameters <- yaml.load_file(param.file)

cape.obj <- new(
  "Cape",
  pheno = pheno,
  chromosome = chromosome,
  marker_num = marker_num,
  marker_location = marker_location,
  geno_names = geno_names,
  ref_allele = ref_allele,
  geno = geno,
  parameters = parameters
)

run.cape(cape.obj, cape.parameter.file, p.or.q = 0.05, n.cores = n.cores,
         run.singlescan = TRUE, run.pairscan = TRUE, error.prop.coef = TRUE,
         error.prop.perm = TRUE, verbose = TRUE)