#' Calculate empirical p-values
#' 
#' This function uses ecdf to calculate empirical p values
#' given a null distribution and an observed distribution
calc.emp.p <- function(obs.dist, null.dist){
  p.fun <- ecdf(null.dist)
  emp.p <- 1-unlist(lapply(obs.dist, p.fun))
  return(emp.p)
}