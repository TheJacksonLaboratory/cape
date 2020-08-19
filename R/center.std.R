#' This function mean centers and standardizes a vector
#'
#' @param v The vector to be mean-centered and standardized
#'
#' @return A mean-centered and standardized vector
#'
#' @export



center.std <- function(v){
  mean.v <- mean(v, na.rm = TRUE)
  centered <- v - mean.v
  sd.v <- sd(v, na.rm = TRUE)
  final.v <- centered/sd.v
  return(final.v)
}