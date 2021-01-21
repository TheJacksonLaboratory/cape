#' Rank Z normalize
#' 
#' Rank Z normalize data in a vector. We rank the data in 
#' x, divide by n-1, and return the quantiles of the normal
#' distribution using qnorm. This transformation is also 
#' called an inverse normal transformation.
#' 
#' @param x data
#' @param jitter boolean, default is FALSE.  
#' 
#' @return Returns a numeric vector with a normal distribution
#'
#' @importFrom stats qnorm
#'
#' @keywords internal
#' 
rz_transform <- function (x, jitter = FALSE){
  x = rank(x, na.last = "keep", ties.method = "average") / (length(x) + 1)
  return(qnorm(x))
}
