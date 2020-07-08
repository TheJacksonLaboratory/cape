if(!require(here)){install.packages("here")}

test.data.path <- here("tests/testthat/testdata")
# file.name <- file.path(test.data.path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")
param.file <- file.path(test.data.path, "cape.parameters.chesler.DO.exp.group.yml")

pheno.obj<-readRDS(file.path(test.data.path, "CheslerDO.pheno.RData"))
geno.obj<-readRDS(file.path(test.data.path, "CheslerDO.geno.RData"))

# set all the dotted list names in the pheno.obj to snake case
names(pheno.obj) <- gsub("[.]", "_", names(pheno.obj))

pheno.obj <- cape::remove.unused.markers(pheno.obj, geno.obj)

final.cross <- run.cape(pheno.obj, geno.obj, results.file = "cross.RData", p.or.q = 0.05,
                        n.cores = 4, run.singlescan = TRUE, run.pairscan = TRUE, error.prop.coef = TRUE,
                        error.prop.perm = TRUE, initialize.only = FALSE, verbose = TRUE, run.parallel = FALSE,
                        param.file = param.file, results.path = here("results"))
