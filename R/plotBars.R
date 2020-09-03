#' Plot phenotypic effect for two markers as a bar plot
#' 
#' This internal function is called by 
#' \code{\link{plot.effects}} to generate a 
#' bar plot showing mean trait values for all
#' combinations of genotypes for two markers.
#' This function also indicates the predicted
#' additive effects of the two markers with 
#' a dashed line, as well as the error of the
#' predicted additive effect with an orange box.
#' The true effect of the marker combination is
#' shown with a gray bar. 
#' The effects can centered on the reference
#' allele to better show relative effects of 
#' each allele and allele combination.
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
#' @param ref.centered Whether to center the effects
#' on the reference genotype.
#' @return None

plotBars <- function(phenoV, marker1.vals, marker2.vals, pheno.name, marker1.label,
marker2.label, ymin = NULL, ymax = NULL, error.bars, ref.centered){

		error.bar.width = 0.1
		addline.width = 0.5
		addline.offset = 0.55
		error.bar.lwd = 2
		text.offset = 0.1
		
		if(is.null(marker2.vals)){
			genotypes <- sort(unique(marker1.vals[which(!is.na(marker1.vals))]))
			pheno.vals <- lapply(genotypes, function(x) phenoV[which(marker1.vals == x)])
			pheno.means <- sapply(pheno.vals, function(x) mean(x, na.rm = TRUE))
			if(error.bars == "sd"){
				pheno.error <- sapply(pheno.vals, function(x) sd(x, na.rm = TRUE))
				}
			if(error.bars == "se"){
				pheno.error <- sapply(pheno.vals, function(x) sd(x, na.rm = TRUE)/sqrt(length(x)))
				}
			if(error.bars == "none"){
				pheno.error <- rep(0, length(pheno.vals))
				}
			# pheno.error[which(is.na(pheno.error))] <- 0
			if(is.null(ymin)){ymin <- min(c(pheno.means - pheno.error, 0))}
			if(is.null(ymax)){ymax <- max(pheno.means + pheno.error)}
			plot.height = ymax - ymin
				
			a <- barplot(pheno.means, ylim = c(ymin, ymax*1.1))
			segments(a, pheno.means-pheno.error, a, pheno.means+pheno.error, lwd = error.bar.lwd)
			abline(h = 0)
			#lower bar
			segments((a+error.bar.width), (pheno.means-pheno.error), (a-error.bar.width), pheno.means-pheno.error, lwd = error.bar.lwd)	
			#upper bar
			segments(a+error.bar.width, pheno.means+pheno.error, a-error.bar.width, pheno.means+pheno.error, lwd = error.bar.lwd)	
			text(x = genotypes, y = ymin-(plot.height*0.1), labels = genotypes)
			mtext(pheno.name, side = 3, line = 1.5)

			}else{ #instead if there are two markers

		
	#========================================================================
	# internal functions
	#========================================================================
	
	plot.bars.int <- function(pheno.vals, pheno.error, pheno.name, marker1, marker2, geno.n){
			
			
			if(length(pheno.vals) == 4){
				base.val = pheno.vals[1]
				}else{
				base.val = 0	
				}
				
			pred.add <- base.val + ((pheno.vals[2]-base.val) + (pheno.vals[3]-base.val))
			if(is.na(pred.add)){pred.add <- 0}
			pred.error <- pheno.error[2] + pheno.error[3]
	
			if(is.null(ymin)){
				if(error.bars != "none"){
					ymin <- min(c((pheno.vals-pheno.error), pred.add, 0), na.rm = TRUE)
					}else{
					ymin <- min(c(pheno.vals, pred.add, 0), na.rm = TRUE)	
					}
				}
			if(is.null(ymax)){
				if(error.bars != "none"){
					ymax <- max(c((pheno.vals+pheno.error), pred.add), na.rm = TRUE)
					}else{
					ymax <- max(c(pheno.vals, pred.add), na.rm = TRUE)	
					}
				}		
			
			
			full.range = ymax-ymin
			geno.n.y <- ymax*1.1
			ymax <- ymax*1.1
			
			label.y <- min(pheno.vals)-full.range*0.1
			a <- barplot(pheno.vals, ylim = c(ymin, ymax), axes = FALSE, 
			xlim = c(0, length(pheno.vals)*1.5), cex.names = 1.5, 
			main = pheno.name, names.arg = NA)
			axis(2, cex.axis = 1.5)
			abline(h = 0)
			par(xpd = TRUE)
			text.cex = 1.5
			text(x = 0, y = ymin-(full.range*0.1), labels = marker1.label, adj = 1, 
			cex = text.cex)
			text(x = 0, y = ymin-(full.range*0.22), labels = marker2.label, adj = 1, 
			cex = text.cex)
			text(x = c(a), y = ymin-(full.range*0.1), labels = c(min.geno, max.geno, 
			min.geno, max.geno), cex = text.cex)
			text(x = c(a), y = ymin-(full.range*0.22), labels = c(min.geno, min.geno, 
			max.geno, max.geno), cex = text.cex)
			if(error.bars != "none"){
				segments(a, pheno.vals-pheno.error, a, pheno.vals+pheno.error, 
				lwd = error.bar.lwd)
				#lower bar
				segments((a+error.bar.width), (pheno.vals-pheno.error), 
				(a-error.bar.width), pheno.vals-pheno.error, lwd = error.bar.lwd)	
				#upper bar
				segments(a+error.bar.width, pheno.vals+pheno.error, a-error.bar.width,
				pheno.vals+pheno.error, lwd = error.bar.lwd)	
				}
			segments(a[length(a)]+addline.width, pred.add, a[length(a)]-addline.width, 
			pred.add, lty = 2, lwd = 2)
			poly.x <- c(a[length(a)]-addline.width, a[length(a)]+addline.width)
			poly.y <- c(pred.add-pred.error, pred.add+pred.error)
			polygon(x = c(poly.x, rev(poly.x)), y = rep(poly.y, each = 2), 
			col = rgb(253/256,192/256,134/256, alpha = 0.5))
			arrows(a[length(a)]+addline.width+addline.offset, pred.add, 
			a[length(a)]+addline.offset, pred.add, lty = 1, lwd = 2, length = 0.1)
			text(x = a[length(a)]+addline.width+addline.offset+text.offset, y = pred.add, 
			labels = "Additive\nPrediction", adj = 0)
			text(x = a[,1], y = rep(geno.n.y, length(a[,1])), labels = geno.n)
			par(xpd = FALSE)
			mtext(pheno.name, side = 2, line = 2.5)
			}			

		#========================================================================

		
		#the genotype combinations used are based on the
		#maximum genotype present. We want to show the 
		#maximum effect for the "mutant" phenotype.
		min.geno <- min(c(marker1.vals, marker2.vals), na.rm = TRUE)
		max.geno <- max(c(marker1.vals, marker2.vals), na.rm = TRUE)
				
				
		w.w <- intersect(which(marker1.vals == min.geno), which(marker2.vals == min.geno))
		w.m <- intersect(which(marker1.vals == max.geno), which(marker2.vals == min.geno))
		m.w <- intersect(which(marker1.vals == min.geno), which(marker2.vals == max.geno))
		m.m <- intersect(which(marker1.vals == max.geno), which(marker2.vals == max.geno))
				
		pheno.list <- vector(mode = "list", length = 4)
		names(pheno.list) <- c(paste(min.geno, min.geno, sep = "/"), 
				paste(min.geno, max.geno, sep = "/"), 
				paste(max.geno, min.geno, sep = "/"), 
				paste(max.geno, max.geno, sep = "/"))
	
		pheno.list[[1]] <- phenoV[w.w]
		pheno.list[[2]] <- phenoV[w.m]
		pheno.list[[3]] <- phenoV[m.w]
		pheno.list[[4]] <- phenoV[m.m]
											
		geno.n <- sapply(pheno.list, length)
		pheno.means <- sapply(pheno.list, function(x) mean(x, na.rm = TRUE))

		if(error.bars == "se"){
			pheno.error <- sapply(pheno.list, function(x) sd(x, na.rm = TRUE)/sqrt(length(x)))
			}
		if(error.bars == "sd"){
			pheno.error <- sapply(pheno.list, function(x) sd(x, na.rm = TRUE))
			}
		if(error.bars == "none"){
			pheno.error <- rep(0, length(pheno.list))
		}


		if(!ref.centered){
			plot.bars.int(pheno.vals = pheno.means, pheno.error = pheno.error, pheno.name = pheno.name, marker1 = marker1.vals, marker2 = marker2.vals, geno.n = geno.n)
			}else{
			#also make a plot using B6 as the reference point
			#center on B6 genotype
			ref.pheno <- pheno.means-pheno.means[1]
			plot.bars.int(pheno.vals = ref.pheno, pheno.error, pheno.name = pheno.name, 
			marker1 = marker1, marker2 = marker2, geno.n = geno.n)
			}		
		}
	} #end case for two markers