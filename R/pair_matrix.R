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
#' @param self_pairs A logical value indicating whether self pairs
#' should be included. If TRUE, a,a will be returned. Otherwise it
#' will be excluded.
#' 
#' @return A two-column matrix containing pairs of entries in elements.
#' @keywords internal
#'
pair_matrix <- function(elements, ordered = FALSE, self_pairs = FALSE){
  
  num_elements <- length(elements)
  
  x_mat <- matrix(elements, ncol = num_elements, nrow = num_elements, byrow = TRUE)
  y_mat <- matrix(elements, ncol = num_elements, nrow = num_elements, byrow = FALSE)
  
  if(ordered){
    
    if(self_pairs){
      upper_x <- c(x_mat[upper.tri(x_mat, diag = TRUE)], x_mat[lower.tri(x_mat, diag = FALSE)])
      upper_y <- c(y_mat[upper.tri(y_mat, diag = TRUE)], y_mat[lower.tri(y_mat, diag = FALSE)])
    }else{
      upper_x <- c(x_mat[upper.tri(x_mat, diag = FALSE)], x_mat[lower.tri(x_mat, diag = FALSE)])
      upper_y <- c(y_mat[upper.tri(y_mat, diag = FALSE)], y_mat[lower.tri(y_mat, diag = FALSE)])
    }
    
  }else{
    
    if(self_pairs){
      upper_x <- x_mat[upper.tri(x_mat, diag = TRUE)]
      upper_y <- y_mat[upper.tri(y_mat, diag = TRUE)]
    }else{
      upper_x <- x_mat[upper.tri(x_mat, diag = FALSE)]
      upper_y <- y_mat[upper.tri(y_mat, diag = FALSE)]
    }
  }
  
  
  pairs_mat <- cbind(upper_y, upper_x)
  colnames(pairs_mat) <- NULL
  return(pairs_mat)
  
}
