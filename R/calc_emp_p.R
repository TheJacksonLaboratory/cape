#' Calculate empirical p-values
#' 
#' This function uses ecdf to calculate empirical p values
#' given a null distribution and an observed distribution
#' 
#' @param obs_dist The observed distribution
#' @param null_dist The null distribution 
#'
#' @return An empirical p value for each observed value
#' @importFrom stats ecdf

calc_emp_p <- function(obs_dist, null_dist){
  p_fun <- ecdf(null_dist)
  emp_p <- 1-unlist(lapply(obs_dist, p_fun))
  return(emp_p)
}