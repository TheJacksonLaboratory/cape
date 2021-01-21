#' This function mean centers and standardizes a vector
#'
#' @param v The vector to be mean-centered and standardized
#'
#' @return A mean-centered and standardized vector
#'
#' @keywords internal


center_std <- function(v){
  mean_v <- mean(v, na.rm = TRUE)
  centered <- v - mean_v
  sd_v <- sd(v, na.rm = TRUE)
  final_v <- centered/sd_v
  return(final_v)
}