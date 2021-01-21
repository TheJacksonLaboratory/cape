#' Drop the 3rd dimension of an array using a summary 
#' function, e.g., min max, mean
#' 
#' This function takes in a 3D array of values, effects or t_stats, etc.
#' and returns a 2D matrix in which each entry contains a single value
#' representing all corresponding entries in the 3rd dimension. This
#' number is determined by the function specified by the user. It could be min,
#' max, mean, etc.
#'
#' @param arrayX 3D array of values, effects or t_stats, etc.
#' @param margin1 dimension that will end up as the columns of the final matrix
#' @param margin2 dimension that will end up as the rows of the final matrix
#' @param slice_fun summary function. The function by which to summarize 
#' the remaining dimension e.g., \code{function(x) mean(x, na.rm = TRUE)}
#'
#' @return A 2D matrix summarizing the input array
#' 
#'
#' @keywords internal
#' 
flatten_array <- function(arrayX, margin1, margin2, slice_fun){
	
	array_dim <- unlist(dim(arrayX))
	final_col <- array_dim[margin2]
	
	slice_dim <- array_dim[-margin1]
	new_margin2 <- which(slice_dim == final_col)
	
	apply_fun <- match.fun(slice_fun)

	flattened_mat <- apply(arrayX, margin1, function(slice) apply(slice, new_margin2, apply_fun))

	return(flattened_mat)


}
