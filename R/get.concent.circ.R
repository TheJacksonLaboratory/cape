#' Generate list of concentric circles
#' 
#' This is an internal function used by \code{\link{PlotNetwork}}
#' and calling get.circle to generate concentric circles starting
#' at a specified radius and with a specified gap between them.
#' These circles are used for plotting traits and chromosomes
#' in concentric circle plots.
#'
#' @param list.names A character vector indicating labels for each circle
#' There will be one circle generated for each element in list.names.
#' @param start.rad The radius for the innermost circle.
#' @param The amount of space between each circle.
#' 
#' @return A list of x,y coordinates as generated
#' by the function \code{\link{get.circle}} and named with 
#' list.names


get.concent.circ <- function(list.names, start.rad = 1.05, gap.rad = 0.05){
  trait.circ <- vector(mode = "list", length = length(list.names))
  names(trait.circ) <- list.names
  for(i in 1:length(list.names)){
    effect.rad = get.circle(start.rad)
    trait.circ[[i]] <- effect.rad
    start.rad <- start.rad + gap.rad
  }
  return(trait.circ)
}
