#' Exponential color function
#' 
#' This function maps colors onto an exponential
#' function for use in colors_from_values.
#' 
#' @param x_min The minimum value assigned to a color
#' @param x_max The maximum value assigned to a color
#' @param steepness The steepness of the exponential function
#' @param num_cols The number of colors to generate in the ramp
#' 
#' @return Numeric vector with x_min:x_max remapped onto an 
#' exponential curve.
#' @keywords internal

exp_color_fun <-
function(x_min, x_max, steepness = 1, num_cols = 256){
	
	x <- seq(0,1,length.out = num_cols)
	y <- exp(steepness*x)-1
	
	#map the 0-1 interval onto the interval from x.min to x.max
	x_range <- x_max - x_min
	scaled_y <- (y*x_range)/max(y)
	

	shifted_y <- scaled_y + x_min
	
	return(shifted_y)
	
	
}
