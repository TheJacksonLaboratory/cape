#' Plot phenotypic effects for two markers as a heat map
#' 
#' This internal function is called by 
#' \code{\link{plot_effects}} to generate a 
#' heat map showing the effects of genotype on
#' phenotype. This function fits linear models
#' to the markers and traits. It then uses
#' these models to predict trait values at 
#' different genotype combinations in a 2D 
#' grid. It plots these predicted values in
#' a heat map.
#' 
#' @param phenoV A vector of trait values 
#' @param marker1_vals A vector of genotype values 
#' for marker1
#' @param marker2_vals A vector of genotype values
#' for marker2.
#' @param pheno_name A string indicating the name of
#' the trait being plotted.
#' @param marker1_label A string indicating the name
#' of marker1
#' @param marker2_label A string indicating the name
#' of marker2
#' @param bins1 The number of bins for marker1 over 
#' which to predict values of the trait. This can also 
#' be a vector specifying those bins.
#' @param bins2 The number of bins for marker2 over 
#' which to predict values of the trait. This can also 
#' be a vector specifying those bins.
#'
#' @return None
#' 
#' @importFrom stats as.formula lm predict
#' @keywords internal

plot_int_heat <- function(phenoV, marker1_vals, marker2_vals, pheno_name = NULL, 
                        marker1_label = NULL, marker2_label = NULL, bins1 = 50, bins2 = 50){

	has_dash1 <- grep("-", marker1_label)
	if(length(has_dash1) > 0){
		marker1_label <- gsub("-", "_", marker1_label)
	}
	has_dash2 <- grep("-", marker2_label)
	if(length(has_dash2) > 0){
		marker2_label <- gsub("-", "_", marker2_label)
	}
	

    if(length(bins1) == 1){
	    min_marker1 <- min(signif(marker1_vals, 2), na.rm = TRUE)
		max_marker1 <- max(signif(marker1_vals, 2), na.rm = TRUE)
        marker_grid1 <- segment_region(min_marker1, max_marker1, bins1, alignment = "ends")
    }else{
    	marker_grid1 <- bins1
    }
    
    if(length(bins2) == 1){
		min_marker2 <- min(signif(marker2_vals, 2), na.rm = TRUE)
		max_marker2 <- max(signif(marker2_vals, 2), na.rm = TRUE)
        marker_grid2 <- segment_region(min_marker2, max_marker2, bins2, alignment = "ends")
    }else{
       	marker_grid2 <- bins2
   	}
        
    marker1_bins <- bin_vector(signif(marker1_vals, 2), bins = marker_grid1)
    marker2_bins <- bin_vector(signif(marker2_vals, 2), bins = marker_grid2)
         	
	test_df <- cbind(phenoV, marker1_bins, marker2_bins)
	colnames(test_df) <- c(pheno_name, marker1_label, marker2_label)
	test_df <- data.frame(test_df)
	fmla <- paste(pheno_name, "~", marker1_label, "*", marker2_label)
    model <- lm(as.formula(fmla), data = test_df)
	#summary(model)
	
	predict_grid <- cbind(rep(marker_grid1, length(marker_grid2)), rep(marker_grid2, each = length(marker_grid1)))
	colnames(predict_grid) <- c(marker1_label, marker2_label)
	predict_df <- data.frame(predict_grid)
	pred_data <- predict(model, newdata = predict_df)
	pred_mat <- matrix(pred_data, nrow = length(marker_grid2), ncol = length(marker_grid1), 
	byrow = FALSE)
	image_with_text(pred_mat, col_text_rotation = 0, use_pheatmap_colors = TRUE, 
	show_text = FALSE, ylab = marker2_label, xlab = marker1_label, 
	main = pheno_name)


}