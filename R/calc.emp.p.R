#' Calculate empirical p-values
#' 
#' This function uses ecdf to calculate empirical p values
#' given a null distribution and an observed distribution
#' 
#' @param obs.dist The observed distribution
#' @param null.dist The null distribution 
#'
#' @return An empirical p value for each observed value
#' @export
calc.emp.p <- function(obs.dist, null.dist){
  p.fun <- ecdf(null.dist)
  emp.p <- 1-unlist(lapply(obs.dist, p.fun))
  return(emp.p)
}