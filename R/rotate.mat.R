#' orients a matrix for proper display in a plot
#' 
#' This function reassembles a matrix 
#' so that when we use image to plot, 
#' the matrix appears in the image in 
#' the same orientation in which it
#' is printed to the screen
#' are in rows and the phenotypes in columns
#' 
#' @param mat a matrix
#' 
rotate.mat <- function(mat){
  n.row <- dim(mat)[1] #get the number of rows and columns
  n.col <- dim(mat)[2]
  new.mat <- matrix(NA, nrow = n.col, ncol = n.row)
  #reverse each column i of the original
  #matrix and put it into row i in the
  #new matrix
  for(i in 1:n.col){
    new.mat[i,] <- rev(mat[,i])
  }
  colnames(new.mat) <- rownames(mat)
  rownames(new.mat) <- colnames(mat)
  return(new.mat)
}
