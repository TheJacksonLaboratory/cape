#This function mean centers and standardizes a vector
center.std <- function(v){
  mean.v <- mean(v, na.rm = TRUE)
  centered <- v - mean.v
  sd.v <- sd(v, na.rm = TRUE)
  final.v <- centered/sd.v
  return(final.v)
}