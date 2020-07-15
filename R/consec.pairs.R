#' Generate a matrix of consecutive elements
#'
#' Given a vector, this function generates a two-column
#' matrix in which each row contains an element and its
#' consecutive neighbor.
#'
#' @param elements a vector of elements
#'
#' @return a matrix with two columns. The first column
#' is identical to elements with the last element taken
#' off. And second column holds the next consecutive 
#' element.
#'
#' @export
consec.pairs <- function(elements){
  
  pair.mat <- matrix(NA, nrow = (length(elements)-1), ncol = 2)
  pair.mat[,1] <- elements[1:(length(elements)-1)]
  pair.mat[,2] <- elements[2:length(elements)]
  return(pair.mat)
}