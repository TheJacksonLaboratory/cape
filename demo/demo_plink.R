setseed(1234)

if(!require(cape)){install.packages("cape")}

results.path <- here("demo", "demo_PLINK")
data.path <- here("tests", "testthat", "testdata", "demo_PLINK_data")

ped <- file.path(data.path, "test.ped")
map <- file.path(data.path, "test.map")
pheno <- file.path(data.path, "test.pheno")
out <- file.path(data.path, "test.csv")
param.file <- file.path(results.path, "plink.parameters.yml")

#unlink(out) #remove the out file if it exists
cross.obj <- plink2cape(ped, map, pheno, out)

data.obj <- cross.obj$data.obj
geno.obj <- cross.obj$geno.obj$geno

final.cross <- run.cape(pheno.obj = data.obj, geno.obj, results.file = "plink.RData", 
	p.or.q = 0.05, verbose = TRUE, param.file = param.file, 
	results.path = results.path)
