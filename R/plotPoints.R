#' Plot phenotypic effect for two markers as points
#' 
#' This internal function is called by 
#' \link{\code{plot.effects}} to generate a 
#' plot showing trait the effects of genotype
#' on phenotype. If marker2.vals is NULL, this
#' plot shows the main effect of marker1. 
#' Otherwise it shows effects of the two alleles
#' together. 
#' In this plot all individual values are plotted
#' as points using color to separate the genotype
#' combinations. The means for each group are indicated
#' by line segments.
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
#'
#' @return None

plotPoints <- function(phenoV, marker1.vals, marker2.vals, pheno.name, marker1.label, marker2.label, ymin = NULL, ymax = NULL){
	
	
	mean.bar.width = 0.15
	jitter.factor = 0.1
	upper.plot.buffer = 0.5
	lower.plot.buffer = 0.2


	if(is.null(ymin)){
		ymin <- min(phenoV, na.rm = TRUE)	
		}
	if(is.null(ymax)){
		ymax <- max(phenoV, na.rm = TRUE)
		}

	if(is.null(marker2.vals)){

		#get the genotype values
		genotypes <- sort(unique(marker1.vals[which(!is.na(marker1.vals))]))
		geno.col <- colors.from.values(1:length(genotypes), use.pheatmap.colors = TRUE)
		#barplot(rep(1, length(geno.col)), col = geno.col)
		
		stripchart(phenoV~as.factor(marker1.vals), vertical = TRUE, method = "jitter", 
		jitter = jitter.factor, pch = 1, col = geno.col[1:length(genotypes)], 
		xlab = marker1.label, ylab = pheno.name, ylim = c(ymin, ymax), main = pheno.name)
		
		xlim = c(min(as.numeric(genotypes))-lower.plot.buffer, max(as.numeric(genotypes))+upper.plot.buffer)
		
		mean.test <- boxplot(phenoV~as.factor(marker1.vals), plot = FALSE)
		
		segments(x0 = (1:length(genotypes) - mean.bar.width), y0 = mean.test[[1]][3,], x1 = (1:length(genotypes) + mean.bar.width), col = "black", lwd = 3)
		}
	
	
	if(!is.null(marker2.vals)){

		#get the genotype values for marker 2
		genotypes <- sort(unique(marker2.vals[which(!is.na(marker2.vals))]))
		geno.col <- colors.from.values(1:length(genotypes), use.pheatmap.colors = TRUE)

		ind.cols <- rep(NA, length(phenoV))
		for(g in 1:length(genotypes)){
			ind.cols[which(marker2.vals == genotypes[g])] <- geno.col[g]
			}

		plot(jitter(marker1.vals, factor = jitter.factor*5),
		phenoV, col = ind.cols, xlim = c(min(as.numeric(genotypes))-lower.plot.buffer,
		max(as.numeric(genotypes))+upper.plot.buffer), axes = FALSE, xlab = marker1.label,
		ylab = pheno.name, pch = 1, ylim = c(ymin, ymax), main = pheno.name)
		axis(1, at = as.numeric(genotypes), labels = genotypes)
		axis(2)
		for(g in 1:length(genotypes)){
			errors <- get.interaction.error(marker1.vals, marker2.vals, phenoV, 
			error.type = "se")
			segments(x0 = (as.numeric(genotypes) - mean.bar.width), y0 = errors$means[g,], x1 = (as.numeric(genotypes) + mean.bar.width), col = geno.col[g], lwd = 3)
			}
		
		par(xpd = TRUE)
		legend(x = max(marker1.vals, na.rm = TRUE)*1.15, 
		y = max(phenoV, na.rm = TRUE), 
		legend = genotypes, col = geno.col, pch = 1, title = marker2.label)
		par(xpd = FALSE)

		
		} #end case for if there are two markers
		
	} #end function