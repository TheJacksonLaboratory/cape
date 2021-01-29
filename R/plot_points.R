#' Plot phenotypic effect for two markers as points
#' 
#' This internal function is called by 
#' \code{\link{plot_effects}} to generate a 
#' plot showing trait the effects of genotype
#' on phenotype. If marker2_vals is NULL, this
#' plot shows the main effect of marker1. 
#' Otherwise it shows effects of the two alleles
#' together. 
#' In this plot all individual values are plotted
#' as points using color to separate the genotype
#' combinations. The means for each group are indicated
#' by line segments.
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
#'
#' @importFrom graphics stripchart boxplot segments
#'
#' @return None
#' @keywords internal

plot_points <- function(phenoV, marker1_vals, marker2_vals, pheno_name, marker1_label, marker2_label, ymin = NULL, ymax = NULL){
	
	oldPar <- par(no.readonly = TRUE)
	on.exit(oldPar)

	geno_pal <- brewer.pal(8, "Set2")

	mean_bar_width = 0.15
	jitter_factor = 0.1
	upper_plot_buffer = 0.5
	lower_plot_buffer = 0.2


	if(is.null(ymin)){
		ymin <- min(phenoV, na.rm = TRUE)	
	}
	if(is.null(ymax)){
		ymax <- max(phenoV, na.rm = TRUE)
	}

	if(is.null(marker2_vals)){

		#get the genotype values
		genotypes <- sort(unique(marker1_vals[which(!is.na(marker1_vals))]))
		#geno_col <- colors_from_values(1:length(genotypes), use_pheatmap_colors = TRUE)
		geno_col <- geno_pal[1:length(genotypes)]
		#barplot(rep(1, length(geno_col)), col = geno_col)
		
		stripchart(phenoV~as.factor(marker1_vals), vertical = TRUE, method = "jitter", 
		jitter = jitter_factor, pch = 1, col = geno_col[1:length(genotypes)], 
		xlab = marker1_label, ylab = pheno_name, ylim = c(ymin, ymax), main = pheno_name)
		
		xlim = c(min(as.numeric(genotypes))-lower_plot_buffer, max(as.numeric(genotypes))+upper_plot_buffer)
		
		mean_test <- boxplot(phenoV~as.factor(marker1_vals), plot = FALSE)
		
		segments(x0 = (1:length(genotypes) - mean_bar_width), y0 = mean_test[[1]][3,], x1 = (1:length(genotypes) + mean_bar_width), col = "black", lwd = 3)
	}
	
	
	if(!is.null(marker2_vals)){

		#get the genotype values for marker 2
		genotypes <- sort(unique(marker2_vals[which(!is.na(marker2_vals))]))
		#geno_col <- colors_from_values(1:length(genotypes), use_pheatmap_colors = TRUE)
		geno_col <- geno_pal[1:length(genotypes)]

		ind_cols <- rep(NA, length(phenoV))
		for(g in 1:length(genotypes)){
			ind_cols[which(marker2_vals == genotypes[g])] <- geno_col[g]
		}

		plot(jitter(marker1_vals, factor = jitter_factor*5),
			phenoV, col = ind_cols, xlim = c(min(as.numeric(genotypes))-lower_plot_buffer,
			max(as.numeric(genotypes))+upper_plot_buffer), axes = FALSE, xlab = marker1_label,
			ylab = pheno_name, pch = 1, ylim = c(ymin, ymax), main = pheno_name)
		axis(1, at = as.numeric(genotypes), labels = genotypes)
		axis(2)
		for(g in 1:length(genotypes)){
			errors <- get_interaction_error(marker1_vals, marker2_vals, phenoV, 
			error_type = "se")
			segments(x0 = (as.numeric(genotypes) - mean_bar_width), y0 = errors$means[g,], x1 = (as.numeric(genotypes) + mean_bar_width), col = geno_col[g], lwd = 3)
		}
		
		par(xpd = TRUE)
		legend(x = max(marker1_vals, na.rm = TRUE)*1.15, 
		y = max(phenoV, na.rm = TRUE), 
		legend = genotypes, col = geno_col, pch = 1, title = marker2_label)
		par(xpd = FALSE)

	} #end case for if there are two markers
		
} #end function