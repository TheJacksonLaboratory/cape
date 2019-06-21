#' This function returns points on a circle given the radius, and center coordinates
get.circle <- function(radius, center.x = 1, center.y = 1, dens = 0.0005){
  t <- seq(0,2*pi, dens)	
  x <- radius*cos(t)+center.x
  y <- radius*sin(t)+center.y
  result <- list(x,y); names(result) <- c("x", "y")
  return(result)
}
