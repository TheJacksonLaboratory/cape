#' get all 2-combinations of elements in a vector
#' 
#' This script makes a matrix with two columns 
#' listing all pairwise combinations of the elements
#' of the input vector.
#' the input vector must have all elements you want
#' pairs of. It can be numeric or character strings
#' For example, to get all pairs of numbers between 
#' 1 and 5, markers = c(1,2,3,4,5)
#' 
#' @param elements a vector containing elements from which to draw all possible pairs
#' @param ordered boolean, TRUE if the elements are ordered
#' @param self.pairs boolean, TRUE if we want to include self-pairs
#' 
pair.matrix <- function(elements, ordered = FALSE, self.pairs = FALSE){
  
  num.elements <- length(elements)
  
  x.mat <- matrix(elements, ncol = num.elements, nrow = num.elements, byrow = TRUE)
  y.mat <- matrix(elements, ncol = num.elements, nrow = num.elements, byrow = FALSE)
  
  if(ordered){
    
    if(self.pairs){
      upper.x <- c(x.mat[upper.tri(x.mat, diag = TRUE)], x.mat[lower.tri(x.mat, diag = FALSE)])
      upper.y <- c(y.mat[upper.tri(y.mat, diag = TRUE)], y.mat[lower.tri(y.mat, diag = FALSE)])
    }else{
      upper.x <- c(x.mat[upper.tri(x.mat, diag = FALSE)], x.mat[lower.tri(x.mat, diag = FALSE)])
      upper.y <- c(y.mat[upper.tri(y.mat, diag = FALSE)], y.mat[lower.tri(y.mat, diag = FALSE)])
    }
    
  }else{
    
    if(self.pairs){
      upper.x <- x.mat[upper.tri(x.mat, diag = TRUE)]
      upper.y <- y.mat[upper.tri(y.mat, diag = TRUE)]
    }else{
      upper.x <- x.mat[upper.tri(x.mat, diag = FALSE)]
      upper.y <- y.mat[upper.tri(y.mat, diag = FALSE)]
    }
  }
  
  
  pairs.mat <- cbind(upper.y, upper.x)
  colnames(pairs.mat) <- NULL
  return(pairs.mat)
  
}


