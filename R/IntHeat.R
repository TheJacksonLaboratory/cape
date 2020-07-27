#' Plot phenotypic effects for two markers as a heat map
#' 
#' This internal function is called by 
#' \link{\code{plot.effects}} to generate a 
#' heat map showing the effects of genotype on
#' phenotype. This function fits linear models
#' to the markers and traits. It then uses
#' these models to predict trait values at 
#' different genotype combinations in a 2D 
#' grid. It plots these predicted values in
#' a heat map.
#' 
#' @param phenoV A vector of trait values 
#' @param marker1.vals A vector of genotype values 
#' for marker1
#' @param marker2.vals A vector of genotype values
#' for marker2.
#' @param pheno.name A string indicating the name of
#' the trait being plotted.
#' @param marker1.label A string indicating the name
#' of marker1
#' @param marker2.label A string indicating the name
#' of marker2
#' @param ymin A numeric value indicating the minimum 
#' y value for the plot. If NULL, it will be calculated
#' based on phenoV and error bars.
#' @param ymax A numeric value indicating the maximum
#' y value for the plot. If NULL, it will be calculated
#' based on phenoV and error bars.
#' @param bins1 The number of bins for marker1 over 
#' which to predict values of the trait. This can also 
#' be a vector specifying those bins.
#' @param bins2 The number of bins for marker2 over 
#' which to predict values of the trait. This can also 
#' be a vector specifying those bins.
#'
#' @return None

IntHeat <- function(phenoV, marker1.vals, marker2.vals, pheno.name = NULL, marker1.label = NULL, marker2.label = NULL, bins1 = 50, bins2 = 50){

    if(length(bins1) == 1){
	    	min.marker1 <- min(signif(marker1.vals, 2), na.rm = TRUE)
		max.marker1 <- max(signif(marker1.vals, 2), na.rm = TRUE)
        marker.grid1 <- segment.region(min.marker1, max.marker1, bins1, alignment = "ends")
    		}else{
    		marker.grid1 <- bins1
    	}
    		    	
    if(length(bins2) == 1){
	    	min.marker2 <- min(signif(marker2.vals, 2), na.rm = TRUE)
		max.marker2 <- max(signif(marker2.vals, 2), na.rm = TRUE)
        marker.grid2 <- segment.region(min.marker2, max.marker2, bins2, alignment = "ends")
    	   }else{
    	    	marker.grid2 <- bins2
   	}
        
    marker1.bins <- bin.vector(signif(marker1.vals, 2), bins = marker.grid1)
    marker2.bins <- bin.vector(signif(marker2.vals, 2), bins = marker.grid2)
         	
	test.df <- cbind(phenoV, marker1.bins, marker2.bins)
	colnames(test.df) <- c(pheno.name, marker1.label, marker2.label)
	test.df <- data.frame(test.df)
	fmla <- paste(pheno.name, "~", marker1.label, "*", marker2.label)
    model <- lm(as.formula(fmla), data = test.df)
	#summary(model)
	
	predict.grid <- cbind(rep(marker.grid1, length(marker.grid2)), rep(marker.grid2, each = length(marker.grid1)))
	colnames(predict.grid) <- c(marker1.label, marker2.label)
	predict.df <- data.frame(predict.grid)
	pred.data <- predict(model, newdata = predict.df)
	pred.mat <- matrix(pred.data, nrow = length(marker.grid2), ncol = length(marker.grid1), byrow = FALSE)
	imageWithText(pred.mat, col.text.rotation = 0, use.pheatmap.colors = TRUE, 
	show.text = FALSE, ylab = marker2.label, xlab = marker1.label, 
	main = pheno.name)


}