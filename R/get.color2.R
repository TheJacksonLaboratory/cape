#' Generate color ramp
#' 
#' This function generates color ramps based on an input
#' color. The colors will range from two shades lighter
#' to two shades darker than the original color, based 
#' on the color gap specified
#' 
#' @param col.name A string indicating a color for the ramp
#' @param col.gap A number indicating the size of the gap
#' between colors
#' @param test.plot Whether to plot colors returned by the function
#' 
#' @return A vector of colors forming a color ramp based on the 
#' input name and numeric gap.

get.color2 <- function(col.name, col.gap = 10, test.plot = FALSE){
	
	col.code <- col2rgb(col.name)	
	col.mat <- matrix(col.code, ncol = 5, nrow = 3, byrow = FALSE)
	shade.seq <- seq((col.gap*-2), (col.gap*2), col.gap)

	for(i in 1:ncol(col.mat)){
		col.mat[,i] <- col.mat[,i] + shade.seq[i]
		}

	#adjust any colors that went below 0 or above 256
	below.col <- which(col.mat < 0, arr.ind = TRUE)[,2]
	if(length(below.col) > 0){
		for(i in 1:length(below.col)){
			min.col <- min(col.mat[,below.col[i]])
			col.diff <- 0 - min.col
			col.mat[,below.col[i]] <- col.mat[,below.col[i]] + col.diff
			}
		}
	
	above.col <- which(col.mat > 256, arr.ind = TRUE)[,2]
	if(length(above.col) > 0){
		for(i in 1:length(above.col)){
			max.col <- max(col.mat[,above.col[i]])
			col.diff <- 256 - max.col
			col.mat[,above.col[i]] <- col.mat[,above.col[i]] + col.diff
			}
		}
	
	
	cols <- apply(col.mat, 2, function(x) rgb(x[1]/256, x[2]/256, x[3]/256))
	
	if(test.plot){
		barplot(rep(1, length(cols)), col = cols)
		}
	
	return(cols)
	}