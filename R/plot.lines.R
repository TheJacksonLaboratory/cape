#' Plot interaction plot for traits and genetic markers
#' 
#' This internal function is called by plot.effects
#' to generate an interaction plot for two markers
#' relating to a trait. If marker2.vals is NULL,
#' the function instead shows the main effect for
#' marker1.
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
#' @param error.bars A string indicating the type of error
#' bars to draw. Can be "sd" for standard deviation, "se"
#' for standard error, or "none".
#'
#' @return None

plot.lines <- function(phenoV, marker1.vals, marker2.vals = NULL, pheno.name, marker1.label, marker2.label, ymin = NULL, ymax = NULL, error.bars){
	
	error.bar.width = 0.15


	#show main effects if marker2 is NULL
	if(is.null(marker2.vals)){
		
		#bin the phenotype into the genotype values, which have already been binned
		geno.bins <- sort(unique(marker1.vals[which(!is.na(marker1.vals))]))
		pheno.bins <- lapply(geno.bins, function(x) phenoV[which(marker1.vals == x)])
		names(pheno.bins) <- geno.bins

		#calculate mean and standard error for each group
		pheno.means <- sapply(pheno.bins, function(x) mean(x, na.rm = TRUE))
		
		if(error.bars == "sd"){
			pheno.error <- sapply(pheno.bins, function(x) sd(x, na.rm = TRUE))
			}
		if(error.bars == "se"){
			pheno.error <- sapply(pheno.bins, function(x) sd(x, na.rm = TRUE)/sqrt(length(x)))
			}
		if(error.bars == "none"){
			pheno.error <- rep(0, length(pheno.bins))
			}

		if(is.null(ymin)){ymin <- min(pheno.means-pheno.error, na.rm = TRUE)}
		if(is.null(ymax)){ymax <- max(pheno.means+pheno.error, na.rm = TRUE)}


		plot(pheno.means, type = "b", ylim = c(ymin, ymax), xlim = c(1,length(geno.bins)),
		axes = FALSE, ylab = "", xlab = marker1.label, main = pheno.name, 
		lwd = 3, cex.lab = 1.5, cex.main = 1.5)
		axis(1, at = 1:length(pheno.means), labels = geno.bins, cex.axis = 2)
		axis(2, cex.axis = 2)
		
		if(error.bars != "none"){
			segments(x0 = 1:length(pheno.means), y0 = pheno.means-pheno.error, y1 = pheno.means+pheno.error, lwd = 2)
			}
		}
		
		#if we have values for two markers
		if(!is.null(marker2.vals)){
			
			marker.geno.bins <- list(sort(unique(marker1.vals)), sort(unique(marker2.vals)))
			errors <- get.interaction.error(marker1.vals, marker2.vals, phenoV, 
			error.type = error.bars)	
			
					
			if(is.null(ymin)){
				all.errors <- errors[[2]]
				all.errors[which(is.na(all.errors))] <- 0
				ymin <- min(errors[[1]]-all.errors, na.rm = TRUE)
				}
			if(is.null(ymax)){
				all.errors <- errors[[2]]
				all.errors[which(is.na(all.errors))] <- 0
				ymax <- max(errors[[1]]+all.errors, na.rm = TRUE)
				}

			ylim <- c(ymin, ymax)
			
			not.na <- which(!is.na(phenoV))
			
			interaction.plot(marker1.vals[not.na], marker2.vals[not.na], phenoV[not.na], 
			xlab = marker1.label, trace.label = marker2.label, ylab = "", lwd = 3, 
			cex.lab = 1.5, main = pheno.name, axes = FALSE, ylim = ylim, 
			fixed = TRUE, cex.main = 1.5)

			plot.lim <- par("usr")
			min.area = plot.lim[1]*1.17; max.area = plot.lim[2]*0.85
			axis(1, at = c(seq(min.area, max.area, (max.area-min.area)/(length(marker.geno.bins[[2]])-1))), labels = marker.geno.bins[[2]], cex.axis = 2); axis(2, cex.axis = 2)
			
			if(error.bars != "none"){
				for(i in 1:length(errors$means[,1])){
					segments(1:length(colnames(errors$means)), (errors$means[i,]+errors[[2]][i,]), 1:length(colnames(errors$means)), (errors$means[i,]-errors[[2]][i,]))
					}
				}
		}#end case for if there are two markers
		
	}#end function
		