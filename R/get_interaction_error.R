#' Get error bars for interaction plot
#' 
#' This function gets error bars for an interaction plot
#' it takes in the x, y, and trace factor of the interaction
#' plot and returns a vector indicating how much should
#' default error type is standard error
#' 
#' @param x The vector whose levels will form the x axis
#' @param y The vector whose levels will form the y axis, 
#' the response vector i.e.
#' @param trace Another vector whose levels will form the traces
#' @param error_type Either "sd" for standard deviation, or
#' "se" for standard error.
#'
#' @return returns a list with two elements. The first element
#' called "means" contains the mean y values for the combinations 
#' of the x and trace variables. The second element, called "se"
#' or "sd" depending on error_type contains the errors for the 
#' same combinations.
#' @keywords internal

get_interaction_error <- function(x, y, trace, error_type = c("sd", "se")){
	
	if(length(grep("e", error_type)) > 0){
		error_type = "se"
		}else{
		error_type = "sd"
		}
	
	x_levels <- levels(as.factor(x))
	y_levels <- levels(as.factor(y))
	
	mean_mat <- matrix(NA, ncol = length(x_levels), nrow = length(y_levels))
	error_mat <- matrix(NA, ncol = length(x_levels), nrow = length(y_levels))
	colnames(mean_mat) <- colnames(error_mat) <- x_levels
	rownames(mean_mat) <- rownames(error_mat) <- y_levels
	
	for(i in 1:length(x_levels)){
		for(j in 1:length(y_levels)){
			group <- intersect(which(x == x_levels[i]), which(y == y_levels[j]))
			if(length(group) > 0){
				mean_mat[j,i] <- mean(trace[group], na.rm = TRUE)
				if(error_type == "se"){
					error_mat[j,i] <- sd(trace[group], na.rm = TRUE)/sqrt(length(trace[group][!is.na(trace[group])]))
					}
				if(error_type == "sd"){
					error_mat[j,i] <- sd(trace[group], na.rm = TRUE)
					}
				}
			}		
		}
	
	results <- list(mean_mat, error_mat)
	names(results) <- c("means", error_type)
	return(results)
	
}