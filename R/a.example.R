# where are the data?
plink.path <- "/Users/emersj/projects/cape/data/plink"

# which file?
file.base <- "test*"

# find the `map` and `ped` files in the globbed directory
map <- NULL
ped <- NULL

for (file.name in Sys.glob(file.path(plink.path, file.base))) {
  file.name <- tolower(file.name)
  if ("ped" == file_ext(file.name)) {
    ped <- file.name
  } else if ("map" == file_ext(file.name)) {
    bed <- file.name
  }
}

# drop out if the files are missing
if (is.null(bed) | is.null(ped)) {
  stop("Input path is missing the BED or PED file.")
}

browser()

cape.obj <- new(
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


# run.cape(cape.obj, cape.parameter.file, p.or.q = 0.05, n.cores = n.cores, 
#          run.singlescan = TRUE, run.pairscan = TRUE, error.prop.coef = TRUE, 
#          error.prop.perm = TRUE, verbose = TRUE)