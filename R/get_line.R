#' Get line coordinates
#' 
#' This function generates points along a line
#' whose coordinates are user-defined.
#'
#' @param x0 The x coordinate of the starting point
#' for the line
#' @param y0 The y coordinate of the starting point
#' for the line
#' @param x1 The x coordinate of the ending point
#' for the line
#' @param y1 The y coordinate of the ending point
#' for the line
#' @param dens A numerical value indicating the density
#' of points to define along the line. Small values lead
#' to more densely calculated points.
#' 
#' @return A list with two elements, x and y. These 
#' elements hold the x and y coordinates respectively 
#' for drawing the line.
#' @keywords internal


get_line <- function(x0, y0, x1, y1, dens = 0.0005){
  x_pts <- segment_region(x0, x1, 1/dens, alignment = "ends")
  slope <- (y1-y0)/(x1-x0)
  y_pts <- slope*(x_pts - x0) + y0
  result <- list(x_pts, y_pts); names(result) <- c("x", "y")
  return(result)
}

