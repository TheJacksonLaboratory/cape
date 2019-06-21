#' This function makes a named list of concentric circles starting 
#' at a specified point and using a specified gap
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
