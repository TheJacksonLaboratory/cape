#' Generate coordinates for a circle
#'
#' Given x,y coordinates and a radius, this function
#' generates points on the circumference of the circle.
#' This function is used by PlotNetwork to plot cape
#' results in a circular layout.
#'
#' @param radius A numeric value giving the radius of the
#' circle
#' @param center.x The x coordinate for the center of the circle
#' @param center.y The y coordinate for the center of the circle
#' @param dens A numeric value controlling how many points are returned.
#' Smaller values will return more points along the circumference of the
#' circle
#'
#' @return This function returns a list with two elements x and y.
#' These are the x and y coordinates of the circle. Plotted against
#' each other they will plot a circle.
#' 

get.circle <- function(radius, center.x = 1, center.y = 1, dens = 0.0005){
  t <- seq(0,2*pi, dens)	
  x <- radius*cos(t)+center.x
  y <- radius*sin(t)+center.y
  result <- list(x,y); names(result) <- c("x", "y")
  return(result)
}
