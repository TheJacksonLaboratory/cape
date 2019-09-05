# install.packages("usethis", repos="http://cran.r-project.org")

library(devtools)
library(here)
library(utils)

dest_path = here("../cape_pkg")

build_path <- devtools::build(pkg = here(), path = dest_path, binary = FALSE, quiet = TRUE)
doc <- devtools::document(pkg = here(), roclets=c('rd', 'collate', 'namespace'))
utils::install.packages(build_path, type="source", lib = dest_path, )
