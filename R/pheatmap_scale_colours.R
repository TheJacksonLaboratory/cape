#' pheatmap scale colours found at this link
#' https://cran.r-project.org/package=pheatmap
#' 
#' This function is internal to 'pheatmap' and not exported, but 
#' our code depends on it. We added it here till it becomes exported.
#' 
#' @param mat matrix
#' @param col color palette default set to rainbow(10)
#' @param breaks number of breaks 
#' @param na_col color for na values
#' 
#' @importFrom grDevices rainbow
#' @keywords internal
#' 
pheatmap_scale_colours <- function (mat, col = rainbow(10), breaks = NA, na_col) 
{
  
  scale_vec_colours <- function (x, col = rainbow(10), breaks = NA, na_col) 
  {
    res <- col[as.numeric(cut(x, breaks = breaks, include.lowest = TRUE))]
    res[is.na(res)] <- na_col
    return(res)
  }
  
  mat = as.matrix(mat)
  return(matrix(scale_vec_colours(as.vector(mat), col = col, 
                                  breaks = breaks, na_col = na_col), nrow(mat), ncol(mat), 
                dimnames = list(rownames(mat), colnames(mat))))
}