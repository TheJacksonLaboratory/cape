#' Orients a matrix for proper display in a plot
#' 
#' This function reassembles a matrix 
#' so that when we use image to plot, 
#' the matrix appears in the image in 
#' the same orientation in which it
#' is printed to the screen
#' are in rows and the phenotypes in columns
#' 
#' @param mat a matrix
#' @return The matrix "mat" rotated 90 degrees.
#' @keywords internal
#'

rotate_mat <- function(mat){
  n_row <- dim(mat)[1] #get the number of rows and columns
  n_col <- dim(mat)[2]
  new_mat <- matrix(NA, nrow = n_col, ncol = n_row)
  #reverse each column i of the original
  #matrix and put it into row i in the
  #new matrix
  for(i in 1:n_col){
    new_mat[i,] <- rev(mat[,i])
  }
  colnames(new_mat) <- rownames(mat)
  rownames(new_mat) <- colnames(mat)
  return(new_mat)
}
