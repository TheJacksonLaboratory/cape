#' This function mean centers and standardizes a vector
#'
#' @param v The vector to be mean-centered and standardized
#'
#' @return A mean-centered and standardized vector
#'
#' #' @examples 
#' randV <- runif(10000)
#' hist(randV, main = "Histogram of Uniform Distribution")
#' centV <- center_std(randV)
#' hist(centV, main = "Histogram of Distribution After Mean Centering and Standardizing")
#' @export


center_std <- function(v){
  mean_v <- mean(v, na.rm = TRUE)
  centered <- v - mean_v
  sd_v <- sd(v, na.rm = TRUE)
  final_v <- centered/sd_v
  return(final_v)
}