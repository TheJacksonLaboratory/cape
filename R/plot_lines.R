#' Plot interaction plot for traits and genetic markers
#' 
#' This internal function is called by plot.effects
#' to generate an interaction plot for two markers
#' relating to a trait. If marker2_vals is NULL,
#' the function instead shows the main effect for
#' marker1.
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
#' @param ymin A numeric value indicating the minimum 
#' y value for the plot. If NULL, it will be calculated
#' based on phenoV and error bars.
#' @param ymax A numeric value indicating the maximum
#' y value for the plot. If NULL, it will be calculated
#' based on phenoV and error bars.
#' @param error_bars A string indicating the type of error
#' bars to draw. Can be "sd" for standard deviation, "se"
#' for standard error, or "none".
#'
#' @return None
#' @importFrom stats interaction.plot
#' @keywords internal

plot_lines <- function(phenoV, marker1_vals, marker2_vals = NULL, pheno_name, marker1_label, marker2_label, ymin = NULL, ymax = NULL, error_bars){
	
	oldPar <- par(no.readonly = TRUE)
	on.exit(oldPar)

	error_bar_width = 0.15


	#show main effects if marker2 is NULL
	if(is.null(marker2_vals)){
		
		#bin the phenotype into the genotype values, which have already been binned
		geno_bins <- sort(unique(marker1_vals[which(!is.na(marker1_vals))]))
		pheno_bins <- lapply(geno_bins, function(x) phenoV[which(marker1_vals == x)])
		names(pheno_bins) <- geno_bins

		#calculate mean and standard error for each group
		pheno_means <- sapply(pheno_bins, function(x) mean(x, na.rm = TRUE))
		
		if(error_bars == "sd"){
			pheno_error <- sapply(pheno_bins, function(x) sd(x, na.rm = TRUE))
		}
		if(error_bars == "se"){
			pheno_error <- sapply(pheno_bins, function(x) sd(x, na.rm = TRUE)/sqrt(length(x)))
		}
		if(error_bars == "none"){
			pheno_error <- rep(0, length(pheno_bins))
		}

		if(is.null(ymin)){ymin <- min(pheno_means-pheno_error, na.rm = TRUE)}
		if(is.null(ymax)){ymax <- max(pheno_means+pheno_error, na.rm = TRUE)}


		plot(pheno_means, type = "b", ylim = c(ymin, ymax), xlim = c(1,length(geno_bins)),
		axes = FALSE, ylab = "", xlab = marker1_label, main = pheno_name, 
		lwd = 3, cex.lab = 1.5, cex.main = 1.5)
		axis(1, at = 1:length(pheno_means), labels = geno_bins, cex.axis = 2)
		axis(2, cex.axis = 2)
		
		if(error_bars != "none"){
			segments(x0 = 1:length(pheno_means), y0 = pheno_means-pheno_error, y1 = pheno_means+pheno_error, lwd = 2)
		}
	}
		
	#if we have values for two markers
	if(!is.null(marker2_vals)){
			
		marker_geno_bins <- list(sort(unique(marker1_vals)), sort(unique(marker2_vals)))
		errors <- get_interaction_error(marker1_vals, marker2_vals, phenoV, 
		error_type = error_bars)	
			
		if(is.null(ymin)){
			all_errors <- errors[[2]]
			all_errors[which(is.na(all_errors))] <- 0
			ymin <- min(errors[[1]]-all_errors, na.rm = TRUE)
		}
		if(is.null(ymax)){
			all_errors <- errors[[2]]
			all_errors[which(is.na(all_errors))] <- 0
			ymax <- max(errors[[1]]+all_errors, na.rm = TRUE)
		}

		ylim <- c(ymin, ymax)
			
		not_na <- which(!is.na(phenoV))
			
		interaction.plot(marker1_vals[not_na], marker2_vals[not_na], phenoV[not_na], 
			xlab = marker1_label, trace.label = marker2_label, ylab = "", lwd = 3, 
			cex.lab = 1.5, main = pheno_name, axes = FALSE, ylim = ylim, 
			fixed = TRUE, cex.main = 1.5)

		plot_lim <- par("usr")
		min_area = plot_lim[1]*1.17; max_area = plot_lim[2]*0.85
		axis(1, at = c(seq(min_area, max_area, (max_area-min_area)/(length(marker_geno_bins[[2]])-1))), labels = marker_geno_bins[[2]], cex.axis = 2); axis(2, cex.axis = 2)
			
		if(error_bars != "none"){
			for(i in 1:length(errors$means[,1])){
				segments(1:length(colnames(errors$means)), (errors$means[i,]+errors[[2]][i,]), 1:length(colnames(errors$means)), (errors$means[i,]-errors[[2]][i,]))
			}
		}
	}#end case for if there are two markers
		
}#end function
		