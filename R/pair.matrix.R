#' Get all pairs of elements in a vector
#' 
#' This function makes a matrix with two columns 
#' listing all pairwise combinations of the elements
#' of the input vector.
#'
#' @param elements a vector containing elements from 
#' which to draw all possible pairs. Can be numeric or 
#' strings
#' @param ordered A logical value indicating whether the order of the
#' pairs is important. If TRUE both a,b and b,a will be returned.
#' If FALSE, only a,b will be returned.
#' @param self.pairs A logical value indicating whether self pairs
#' should be included. If TRUE, a,a will be returned. Otherwise it
#' will be excluded.
#' 
#' @return A two-column matrix containing pairs of entries in elements.
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


