## ---- install_cape, eval = FALSE----------------------------------------------
#  install.packages("cape")

## ----load_cape, echo = FALSE, warning = FALSE, error = FALSE, message = FALSE----
set.seed(1234)
library(cape)

## ----read_csv_format----------------------------------------------------------
results_path <- here("demo", "demo_qtl")
data_path <- here("tests", "testthat", "testdata", "demo_qtl_data")
data_file <- file.path(data_path, "NON_NZO_Reifsnyder_pgm_CAPE_num.csv")
param_file <- file.path(results_path, "NON_NZO.parameters.yml")

cross <- read_population(data_file)
cross_obj <- cape2mpp(cross)
obesity_cross <- cross_obj$data_obj
obesity_geno <- cross_obj$geno_obj$geno

## ----read_qtl2, eval = FALSE--------------------------------------------------
#  iron_qtl2 <- read_cross2("https://kbroman.org/qtl2/assets/sampledata/iron/iron.yaml")
#  
#  iron_cape <- qtl2_to_cape(cross = iron_qtl2)
#  data_obj <- iron_cape$data_obj
#  geno_obj <- iron_cape$geno_obj

## ----read_plink, eval = FALSE-------------------------------------------------
#  data_path <- here("tests", "testthat", "testdata")
#  ped <- file.path(data_path, "test.ped")
#  map <- file.path(data_path, "test.map")
#  pheno <- file.path(data_path, "test.pheno")
#  out <- file.path(data_path, "test.csv")
#  param_file <- file.path(results_path, "plink.parameters.yml")
#  
#  cross_obj <- plink2cape(ped, map, pheno, out = "out.csv")
#  
#  data_obj <- cross_obj$data_obj
#  geno_obj <- cross_obj$geno_obj$geno

## ----pheno_hist, fig.width = 7, fig.height = 3--------------------------------
hist_pheno(obesity_cross, pheno_which = c("BW_24", "INS_24", "log_GLU_24"))

## ----qnorm_pheno, fig.width = 7, fig.height = 3-------------------------------
qnorm_pheno(obesity_cross, pheno_which = c("BW_24", "INS_24", "log_GLU_24"))

## ----norm_pheno---------------------------------------------------------------
obesity_cross <- norm_pheno(obesity_cross, mean_center = TRUE)

## ----norm_qq, fig.width = 7, fig.height = 3-----------------------------------
qnorm_pheno(obesity_cross, pheno_which = c("BW_24", "INS_24", "log_GLU_24"))

## ----pheno_cor, fig.width = 5, fig.height = 5---------------------------------
plot_pheno_cor(obesity_cross, pheno_which = c("BW_24", "INS_24", "log_GLU_24"),
color_by = "pgm", group_labels = c("Non-obese", "Obese"))

## ----run_cape-----------------------------------------------------------------
final_cross <- run_cape(obesity_cross, obesity_geno, results_file = "NON_NZO.RData", 
p_or_q = 0.05, verbose = FALSE, param_file = param_file, results_path = results_path)

## ----eigentraits, results = "asis", echo = FALSE, out.width = 4, out.height = 4----
svd_file <- here("demo", "demo_qtl", "svd.jpg")
cat(paste0("![](", svd_file, ")\n"))

## ----single_plot, fig.width = 7, fig.height = 5-------------------------------
singlescan_obj <- readRDS(here("demo", "demo_qtl", "NON_NZO_singlescan.RData"))
plot_singlescan(final_cross, singlescan_obj, line_type = "h", lwd = 2, 
covar_label_size = 1)

## ----reparam_fig, results = "asis", echo = FALSE, fig.width = 3, fig.height = 3----
reparam_file <- here("vignettes", "reparam.png")
cat(paste0("![](", reparam_file, ")\n"))

## ----var_inf, results = "asis", echo = FALSE----------------------------------
var_inf_file <- here("demo", "demo_qtl", "variant_influences.jpg")
cat(paste0("![](", var_inf_file, ")\n"))

## ----circ_net, results = "asis", echo = FALSE---------------------------------
circ_net_file <- here("demo", "demo_qtl", "Network_Circular.jpg")
cat(paste0("![](", circ_net_file, ")\n"))

## ----net_vis, results = "asis", echo = FALSE----------------------------------
net_file <- here("demo", "demo_qtl", "Network_View.jpg")
cat(paste0("![](", net_file, ")\n"))

## ----plotrlang::last_error()_main, fig.width = 7, fig.height = 3--------------
plot_effects(data_obj = final_cross, geno_obj = obesity_geno, 
marker1 = "D15Mit72_B", marker1_label = "Chr15", plot_type = "l", 
error_bars = "se")

## ----plot_int, fig.width = 7, fig.height = 3----------------------------------
plot_effects(data_obj = final_cross, geno_obj = obesity_geno, 
marker1 = "D2Mit120_B", marker2 = "D15Mit72_B", marker1_label = "Chr2",
marker2_label = "Chr15", plot_type = "l", error_bars = "se")

## ----int_bar, fig.width = 7, fig.height = 3-----------------------------------
plot_effects(data_obj = final_cross, geno_obj = obesity_geno, 
marker1 = "D2Mit120_B", marker2 = "D15Mit72_B", marker1_label = "Chr2",
marker2_label = "Chr15", plot_type = "b", error_bars = "se")

