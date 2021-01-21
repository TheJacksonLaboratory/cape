#' Generate color ramp
#' 
#' This function generates color ramps based on an input
#' color. The colors will range from two shades lighter
#' to two shades darker than the original color, based 
#' on the color gap specified
#' 
#' @param col_name A string indicating a color for the ramp
#' @param col_gap A number indicating the size of the gap
#' between colors
#' @param test_plot Whether to plot colors returned by the function
#' 
#' @return A vector of colors forming a color ramp based on the 
#' input name and numeric gap.
#' @keywords internal

get_color2 <- function(col_name, col_gap = 10, test_plot = FALSE){
	
	col_code <- col2rgb(col_name)	
	col_mat <- matrix(col_code, ncol = 5, nrow = 3, byrow = FALSE)
	shade_seq <- seq((col_gap*-2), (col_gap*2), col_gap)

	for(i in 1:ncol(col_mat)){
		col_mat[,i] <- col_mat[,i] + shade_seq[i]
		}

	#adjust any colors that went below 0 or above 256
	below_col <- which(col_mat < 0, arr.ind = TRUE)[,2]
	if(length(below_col) > 0){
		for(i in 1:length(below_col)){
			min_col <- min(col_mat[,below_col[i]])
			col_diff <- 0 - min_col
			col_mat[,below_col[i]] <- col_mat[,below_col[i]] + col_diff
			}
		}
	
	above_col <- which(col_mat > 256, arr.ind = TRUE)[,2]
	if(length(above_col) > 0){
		for(i in 1:length(above_col)){
			max_col <- max(col_mat[,above_col[i]])
			col_diff <- 256 - max_col
			col_mat[,above_col[i]] <- col_mat[,above_col[i]] + col_diff
			}
		}
	
	
	cols <- apply(col_mat, 2, function(x) rgb(x[1]/256, x[2]/256, x[3]/256))
	
	if(test_plot){
		barplot(rep(1, length(cols)), col = cols)
		}
	
	return(cols)
	}