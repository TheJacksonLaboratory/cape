# there is a bug in the curl 4.1 package, use 4.0 instead
install.packages("https://cran.r-project.org/src/contrib/Archive/curl/curl_4.0.tar.gz",repo=NULL,type="source")
install.packages(c("testthat", "here", "devtools", "evd","fdrtool","shape","corpcor","doParallel","foreach","caTools","abind","propagate","tidyr","data.table","RcppEigen","RSQLite","qtl","regress","BiocManager","igraph", "RColorBrewer", "pheatmap"),repos = "http://cran.us.r-project.org")

install.packages(c("qtl2", "qtl2convert"), repos="https://rqtl.org/qtl2cran")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")

BiocManager::install(c('biomaRt'), ask = FALSE)
