#' Find column numbers using column names
#'
#' This is an internal function to find the
#' column numbers of phenotypes when names
#' are put in instead of column numbers.
#'
#' @param data_mat a matrix with column names
#' @param col_which a character string indicating which column should be identified
#' If col_which is numeric the number is returned.
#' @param dim_which A number indicating whether the algorithm should look
#' in rows (1) or columns (2) for the name.
#'
#' @return A numeric value indicating the column number corresponding
#' to the name given by col_which. If col_which was a number, col_which
#' is returned.
#' @keywords internal

get_col_num <- function (data_mat, col_which = NULL, dim_which = 2) 
{
  if (is.null(col_which)) {
    col_num <- 1:dim(data_mat)[dim_which]
  }
  else {
    if (is.numeric(col_which[1])) {
      return(col_which)
    }
    else {
      col_mat <- matrix(col_which, ncol = 1)
      find_col <- function(col_name, data_names) {
        col_num <- which(data_names == col_name)
        if(length(col_num) == 0){
          return(0)
        }
        return(col_num)
      }
      
      col_num <- unlist(apply(col_mat, 1, function(x) find_col(x, dimnames(data_mat)[[dim_which]])))
      didnt_find <- which(col_num == 0)
      if(length(didnt_find) > 0) {
        # message("The following headers are not present:")
        # cat(col_which[didnt_find], sep = "\n")
        col_num <- col_num[-didnt_find]
      }
    }
  }
  return(col_num)
}
