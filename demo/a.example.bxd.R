
if(!require(here)){install.packages("here")}
# if(!require(cape)){install.packages("cape")}

test.data.path <- here("tests/testthat/testdata")
file.name <- file.path(test.data.path, "BxD/BxD.yaml")
param.file <- file.path(test.data.path, "cape.parameters.bxd.yml")

qtl2.bxd <- read_cross2(file.name)

cape.object <- qtl2_to_cape(qtl2.bxd)
data.obj <- cape.object$data.obj
geno.obj <- cape.object$geno.obj 

final.cross <- run.cape(data.obj, geno.obj, results.file = "cross.RData", p.or.q = 0.05,
         n.cores = 4, run.singlescan = TRUE, run.pairscan = TRUE, error.prop.coef = TRUE, 
         error.prop.perm = TRUE, initialize.only = TRUE, verbose = TRUE, run.parallel = FALSE,
         param.file = param.file, results.path = here("results"))
