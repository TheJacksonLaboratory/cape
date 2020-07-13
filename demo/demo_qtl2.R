library(cape)

demo.path <- here("demo", "demo_qtl2")
param.file <- file.path(demo.path, "iron.parameters.yml")

#read in example qtl2 data from remote host
iron.qtl2 <- read_cross2("https://kbroman.org/qtl2/assets/sampledata/iron/iron.yaml")

iron.cape <- qtl2_to_cape(cross = iron.qtl2)
data.obj <- iron.cape$data.obj
geno.obj <- iron.cape$geno.obj 

final.cross <- run.cape(pheno.obj = data.obj, geno.obj, results.file = "iron.RData", 
	p.or.q = 0.05, verbose = TRUE, param.file = param.file, 
	results.path = demo.path)
