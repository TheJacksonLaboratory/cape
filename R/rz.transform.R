#' Rank Z noramalize
#' 
#' Rank Z normalize data in a vector. We rank the data in 
#' x, divide by n-1, and return the quantiles of the normal
#' distribution using qnorm. This transformation is also 
#' called an inverse normal transformation.
#' 
#' @return Returns a numeric vector with a normal distribution
#'
#' @export


rz.transform <- function (x, jitter = FALSE){
  x = rank(x, na.last = "keep", ties.method = "average") / (length(x) + 1)
  return(qnorm(x))
}
