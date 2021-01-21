#' Snap continuous values to a grid
#' 
#' This function bins a continuously valued vector
#' based on user-defined bins.
#' It is useful for binning continuously valued genotypes.
#' Each value in the matrix gets shifted to the nearest
#' provided in the argument bins.
#' 
#' @param vectorX A vector of numeric values to bin
#' @param bins A vector of values to snap the values in vectorX to.
#' 
#' @return A vector the same length as vectorX in which 
#' each value in vectorX has been sent to the nearest value
#' in bins. For example, if bins is c(0, 0.5, 1), and vectorX
#' contains a 0.49. That 0.49 value will be sent to 0.5.
#' The returned vector contains only the values specified in bins.
#'
#' @keywords internal

bin_vector <- function(vectorX, bins = seq(0,1,0.5)){
	
	dist_mat <- apply(matrix(bins, ncol = 1), 1, function(x) x - vectorX)	
	binned_vector <- apply(dist_mat, 1, function(x) bins[which(abs(x) == min(abs(x)))[1]])
	return(binned_vector)	
	
}