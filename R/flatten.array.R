#' drop the 3rd dimension of an array using a summary function, e.g., min max, mean
#' 
#' This function takes in a 3D array of values, effects or t.stats, etc.
#' and returns a 2D matrix in wich each entry contains a single value
#' representing all corresponding entries in the 3rd dimension. This
#' number is determined by the function the user uses. It could be min,
#' max, mean, etc.
#'
#' @param arrayX 3D array of values, effects or t.stats, etc.
#' @param margin1 dimension/axis of the rows (usually 1)
#' @param margin2 dimension/axis of the columns (usually 2)
#' @param slice.fun summary function, e.g., \code{function(x) mean(x, na.rm = TRUE)}
#'
#' @return an updated cape object
#'
#' @export
flatten.array <- function(arrayX, margin1, margin2, slice.fun){
	
	array.dim <- unlist(dim(arrayX))
	final.col <- array.dim[margin2]
	
	slice.dim <- array.dim[-margin1]
	new.margin2 <- which(slice.dim == final.col)
	
	apply.fun <- match.fun(slice.fun)

	flattened.mat <- apply(arrayX, margin1, function(slice) apply(slice, new.margin2, apply.fun))

	return(flattened.mat)


}