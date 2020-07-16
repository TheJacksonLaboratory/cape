set.seed(1234)

library(cape)

demo.path <- here("demo", "demo_qtl2")
data.path <- here("tests", "testthat", "testdata", "demo_qtl2_data")

data.file <- file.path(data.path, "iron.yaml")
param.file <- file.path(demo.path, "iron.parameters.yml")

#read in example qtl2 data from remote host
iron.qtl2 <- read_cross2(data.file)

iron.cape <- qtl2_to_cape(cross = iron.qtl2)
data.obj <- iron.cape$data.obj
geno.obj <- iron.cape$geno.obj 

final.cross <- run.cape(pheno.obj = data.obj, geno.obj, results.file = "iron.RData", 
	p.or.q = 0.05, verbose = TRUE, param.file = param.file, 
	results.path = demo.path)
