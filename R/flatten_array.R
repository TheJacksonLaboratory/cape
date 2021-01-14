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
#' @param dimension that will end up as the columns of the final matrix
#' @param dimension that will end up as the rows of the final matrix
#' @param slice_fun summary function. The function by which to summarize 
#' the remaining dimension e.g., \code{function(x) mean(x, na.rm = TRUE)}
#'
#' @return A 2D matrix summarizing the input array
#' 
#' @examples 
#' rand_array <- array(1:12, dim = c(3,4,5))
#' 
#' #flatten by preserving the first and third dimensions
#' #final matrix has dim[3] rows and dim[1] columns
#' #take the mean across the second dimension
#' flat_mat1 <- flatten_array(rand_array, 1, 3, "mean")
#' print(flat_mat1)
#' 
#' #flatten by preserving the first and third dimensions
#' #final matrix has dim[1] rows and dim[3] columns
#' #take the median across the second dimension
#' flat_mat2 <- flatten_array(rand_array, 3, 1, "median")
#' print(flat_mat2)
#' 
#' #flatten by preserving the first and second dimensions
#' #final matrix has dim[2] rows and dim[1] columns
#' #take the maximum across the second dimension
#' flat_mat3 <- flatten_array(rand_array, 1, 2, "max")
#' print(flat_mat3)
#'
#' @export
flatten_array <- function(arrayX, margin1, margin2, slice_fun){
	
	array_dim <- unlist(dim(arrayX))
	final_col <- array_dim[margin2]
	
	slice_dim <- array_dim[-margin1]
	new_margin2 <- which(slice_dim == final_col)
	
	apply_fun <- match.fun(slice_fun)

	flattened_mat <- apply(arrayX, margin1, function(slice) apply(slice, new_margin2, apply_fun))

	return(flattened_mat)


}
