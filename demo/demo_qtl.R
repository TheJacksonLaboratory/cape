setseed(1234)

if(!require(cape)){install.packages("cape")}

results.path <- here("demo", "demo_qtl")
data.path <- here("tests", "testthat", "testdata", "demo_qtl_data")

data.file <- file.path(data.path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")
param.file <- file.path(results.path, "NON_NZO.parameters.yml")

cross <- read.population(data.file)
cross.obj <- cape2mpp(cross)
data.obj <- cross.obj$data.obj
geno.obj <- cross.obj$geno.obj$geno

final.cross <- run.cape(data.obj, geno.obj, results.file = "NON_NZO.RData", 
	p.or.q = 0.05, verbose = TRUE, param.file = param.file, 
	results.path = results.path)
