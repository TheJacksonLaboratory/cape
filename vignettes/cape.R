## ----load_cape, echo = FALSE, warning = FALSE, error = FALSE, message = FALSE----
set.seed(1234)
library(cape)

## ----read_csv_format----------------------------------------------------------
results.path <- here("demo", "demo_qtl")
data.path <- here("tests", "testthat", "testdata", "demo_qtl_data")
data.file <- file.path(data.path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")
param.file <- file.path(results.path, "NON_NZO.parameters.yml")

cross <- read.population(data.file)
cross.obj <- cape2mpp(cross)
obesity.cross <- cross.obj$data.obj
obesity.geno <- cross.obj$geno.obj$geno

## ----read_qtl2, eval = FALSE--------------------------------------------------
#  iron.qtl2 <- read_cross2("https://kbroman.org/qtl2/assets/sampledata/iron/iron.yaml")
#  
#  iron.cape <- qtl2_to_cape(cross = iron.qtl2)
#  data.obj <- iron.cape$data.obj
#  geno.obj <- iron.cape$geno.obj

## ----read_plink, eval = FALSE-------------------------------------------------
#  data.path <- here("tests", "testthat", "testdata")
#  ped <- file.path(data.path, "test.ped")
#  map <- file.path(data.path, "test.map")
#  pheno <- file.path(data.path, "test.pheno")
#  out <- file.path(data.path, "test.csv")
#  param.file <- file.path(results.path, "plink.parameters.yml")
#  
#  cross.obj <- plink2cape(ped, map, pheno, out = "out.csv")
#  
#  data.obj <- cross.obj$data.obj
#  geno.obj <- cross.obj$geno.obj$geno

## ----pheno_hist, fig.width = 7, fig.height = 3--------------------------------
histPheno(obesity.cross, pheno.which = c("BW_24", "INS_24", "log_GLU_24"))

## ----qnormPheno, fig.width = 7, fig.height = 3--------------------------------
qnormPheno(obesity.cross, pheno.which = c("BW_24", "INS_24", "log_GLU_24"))

## ----norm_pheno---------------------------------------------------------------
obesity.cross <- norm.pheno(obesity.cross, mean.center = TRUE)

## ----norm_qq, fig.width = 7, fig.height = 3-----------------------------------
qnormPheno(obesity.cross, pheno.which = c("BW_24", "INS_24", "log_GLU_24"))

## ----pheno_cor, fig.width = 5, fig.height = 5---------------------------------
plotPhenoCor(obesity.cross, pheno.which = c("BW_24", "INS_24", "log_GLU_24"),
color.by = "mom", group.labels = c("Non-obese", "Obese"))

## ----run_cape-----------------------------------------------------------------
final.cross <- run.cape(obesity.cross, obesity.geno, 
  results.file = "NON_NZO.RData", p.or.q = 0.05, verbose = FALSE, 
  param.file = param.file, results.path = results.path)

## ----eigentraits, results = "asis", echo = FALSE, out.width = 4, out.height = 4----
svd.file <- here("demo", "demo_qtl", "svd.jpg")
cat(paste0("![](", svd.file, ")\n"))

## ----single_plot, fig.width = 7, fig.height = 5-------------------------------
singlescan.obj <- readRDS(here("demo", "demo_qtl", "NON_NZO.singlescan.RData"))
plotSinglescan(final.cross, singlescan.obj, line.type = "h", lwd = 2, 
covar.label.size = 1)

## ----reparam_fig, results = "asis", echo = FALSE, fig.width = 3, fig.height = 3----
reparam.file <- here("vignettes", "reparam.png")
cat(paste0("![](", reparam.file, ")\n"))

## ----var_inf, results = "asis", echo = FALSE----------------------------------
var.inf.file <- here("demo", "demo_qtl", "variant.influences.jpg")
cat(paste0("![](", var.inf.file, ")\n"))

## ----circ_net, results = "asis", echo = FALSE---------------------------------
circ.net.file <- here("demo", "demo_qtl", "Network.Circular.jpg")
cat(paste0("![](", circ.net.file, ")\n"))

## ----net_vis, results = "asis", echo = FALSE----------------------------------
net.file <- here("demo", "demo_qtl", "Network.View.jpg")
cat(paste0("![](", net.file, ")\n"))

## ----plot_main, fig.width = 7, fig.height = 3---------------------------------
plotEffects(data.obj = final.cross, geno.obj = obesity.geno, 
marker1 = "D15Mit72_B", marker1.label = "Chr15", plot.type = "l", 
error.bars = "se")

## ----plot_int, fig.width = 7, fig.height = 3----------------------------------
plotEffects(data.obj = final.cross, geno.obj = obesity.geno, 
marker1 = "D2Mit120_B", marker2 = "D15Mit72_B", marker1.label = "Chr2",
marker2.label = "Chr15", plot.type = "l", error.bars = "se")

## ----int_bar, fig.width = 7, fig.height = 3-----------------------------------
plotEffects(data.obj = final.cross, geno.obj = obesity.geno, 
marker1 = "D2Mit120_B", marker2 = "D15Mit72_B", marker1.label = "Chr2",
marker2.label = "Chr15", plot.type = "b", error.bars = "se")

