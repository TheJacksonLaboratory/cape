# #!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# 
# library(devtools)
# library(testthat)
# library(here)
# 
# cat(args[1])
# cat(" ")
# cat(print("cape" %in% (.packages())))
# 
# 
# # set output xml file
# options(testthat.output_file = here("tests/test-out.xml"))
# 
# # run tests
# testthat::test_dir(here("tests/testthat"), env=.GlobalEnv, reporter = "junit")  #, filter = 'run')
