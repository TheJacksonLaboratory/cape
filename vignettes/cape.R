## ---- install_cape, eval = FALSE----------------------------------------------
#  install.packages("cape")

## ----load_cape, echo = FALSE, warning = FALSE, error = FALSE, message = FALSE----
set.seed(1234)
library(cape)

## ----read_csv_format----------------------------------------------------------
results_path <- here::here("demo", "demo_qtl")
data_path <- here::here("tests", "testthat", "testdata", "demo_qtl_data")
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
#  data_path <- here::here("tests", "testthat", "testdata")
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

## ----eigentraits, fig.width = 4, fig.height = 4-------------------------------
plot_svd(final_cross)

## ----single_plot, eval = FALSE------------------------------------------------
#  singlescan_obj <- readRDS(here::here("demo", "demo_qtl", "NON_NZO_singlescan.RData"))
#  plot_singlescan(final_cross, singlescan_obj, line_type = "h", lwd = 2,
#  covar_label_size = 1)

## ----single_plot_fig1, results = "asis", echo = FALSE-------------------------
et1_fig <- here::here("vignettes", "Singlescan_ET1_Standardized.jpg")
cat(paste0("![](", et1_fig, "){width=70%}\n"))

## ----single_plot_fig2, results = "asis", echo = FALSE-------------------------
et2_fig <- here::here("vignettes", "Singlescan_ET2_Standardized.jpg")
cat(paste0("![](", et2_fig, "){width=70%}\n"))

## ----reparam_fig, results = "asis", echo = FALSE------------------------------
reparam_file <- here::here("vignettes", "reparam.png")
cat(paste0("![](", reparam_file, "){width=50%}\n"))

## ----var_inf, fig.height = 6, fig.width = 7-----------------------------------
plot_variant_influences(final_cross, show_alleles = FALSE)

## ----circ_net, fig.height = 6, fig.width = 6----------------------------------
plot_network(final_cross)

## ----net_vis, echo = FALSE, fig.width = 6, fig.height = 6---------------------
plot_full_network(final_cross, zoom = 1.2, node_radius = 0.3, 
    label_nodes = TRUE, label_offset = 0.4, label_cex = 0.5, bg_col = "lightgray", 
    arrow_length = 0.1, layout_matrix = "layout_with_kk", legend_position = "topright", 
    edge_lwd = 1, legend_radius = 2, legend_cex = 0.7, xshift = -1)

## ----plot_effects_line, fig.width = 7, fig.height = 3-------------------------
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

