#' Generate list of concentric circles
#' 
#' This is an internal function used by \code{\link{plot_network}}
#' and calling get_circle to generate concentric circles starting
#' at a specified radius and with a specified gap between them.
#' These circles are used for plotting traits and chromosomes
#' in concentric circle plots.
#'
#' @param list_names A character vector indicating labels for each circle
#' There will be one circle generated for each element in list_names.
#' @param start_rad The radius for the innermost circle.
#' @param gap_rad The amount of space between each circle.
#' 
#' @return A list of x,y coordinates as generated
#' by the function \code{\link{get_circle}} and named with 
#' list_names
#' @keywords internal


get_concent_circ <- function(list_names, start_rad = 1.05, gap_rad = 0.05){
  trait_circ <- vector(mode = "list", length = length(list_names))
  names(trait_circ) <- list_names
  for(i in 1:length(list_names)){
    effect_rad = get_circle(start_rad)
    trait_circ[[i]] <- effect_rad
    start_rad <- start_rad + gap_rad
  }
  return(trait_circ)
}
