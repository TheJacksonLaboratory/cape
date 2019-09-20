install.packages(c("testthat", "here", "devtools", "evd","fdrtool","shape","corpcor","doParallel","foreach","caTools","abind","propagate","tidyr","data.table","RcppEigen","RSQLite","qtl","regress","BiocManager","igraph"),repos = "http://cran.us.r-project.org")

install.packages(c("qtl2", "qtl2convert"), repos="https://rqtl.org/qtl2cran")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")

BiocManager::install(c('biomaRt'), ask = FALSE)