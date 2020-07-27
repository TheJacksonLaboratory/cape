#' Get error bars for interaction plot
#' 
#' This function gets error bars for an interaction plot
#' it takes in the x, y, and trace factor of the interatction
#' plot and returns a vector indicating how much should
#' default error type is standard error
#' 
#' @param x The vector whose levels will form the x axis
#' @param y The vector whose levels will form the y axis, 
#' the response vector i.e.
#' @param trace Another vector whose levels will form the traces
#' @param error.type Either "sd" for standard deviation, or
#' "se" for standard error.
#'
#' @return returns a list with two elements. The first element
#' called "means" contains the mean y values for the combinations 
#' of the x and trace variables. The second element, called "se"
#' or "sd" depending on error.type contains the errors for the 
#' same combinations.

get.interaction.error <- function(x, y, trace, error.type = c("sd", "se")){
	
	if(length(grep("e", error.type)) > 0){
		error.type = "se"
		}else{
		error.type = "sd"
		}
	
	x.levels <- levels(as.factor(x))
	y.levels <- levels(as.factor(y))
	
	mean.mat <- matrix(NA, ncol = length(x.levels), nrow = length(y.levels))
	error.mat <- matrix(NA, ncol = length(x.levels), nrow = length(y.levels))
	colnames(mean.mat) <- colnames(error.mat) <- x.levels
	rownames(mean.mat) <- rownames(error.mat) <- y.levels
	
	for(i in 1:length(x.levels)){
		for(j in 1:length(y.levels)){
			group <- intersect(which(x == x.levels[i]), which(y == y.levels[j]))
			if(length(group) > 0){
				mean.mat[j,i] <- mean(trace[group], na.rm = TRUE)
				if(error.type == "se"){
					error.mat[j,i] <- sd(trace[group], na.rm = TRUE)/sqrt(length(trace[group][!is.na(trace[group])]))
					}
				if(error.type == "sd"){
					error.mat[j,i] <- sd(trace[group], na.rm = TRUE)
					}
				}
			}		
		}
	
	results <- list(mean.mat, error.mat)
	names(results) <- c("means", error.type)
	return(results)
	
}