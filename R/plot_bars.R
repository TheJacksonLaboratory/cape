#' Plot phenotypic effect for two markers as a bar plot
#' 
#' This internal function is called by 
#' \code{\link{plot_effects}} to generate a 
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
#' @param ref_centered Whether to center the effects
#' on the reference genotype.
#' @return None
#' 
#' @importFrom graphics mtext
#' @keywords internal

plot_bars <- function(phenoV, marker1_vals, marker2_vals, pheno_name, marker1_label,
	marker2_label, ymin = NULL, ymax = NULL, error_bars, ref_centered){

		oldPar <- par(no.readonly = TRUE)
		on.exit(oldPar)

		error_bar_width = 0.1
		addline_width = 0.5
		addline_offset = 0.55
		error_bar_lwd = 2
		text_offset = 0.1
		
		if(is.null(marker2_vals)){
			genotypes <- sort(unique(marker1_vals[which(!is.na(marker1_vals))]))
			pheno_vals <- lapply(genotypes, function(x) phenoV[which(marker1_vals == x)])
			pheno_means <- sapply(pheno_vals, function(x) mean(x, na.rm = TRUE))
			if(error_bars == "sd"){
				pheno_error <- sapply(pheno_vals, function(x) sd(x, na.rm = TRUE))
			}
			if(error_bars == "se"){
				pheno_error <- sapply(pheno_vals, function(x) sd(x, na.rm = TRUE)/sqrt(length(x)))
			}
			if(error_bars == "none"){
				pheno_error <- rep(0, length(pheno_vals))
			}
			# pheno_error[which(is.na(pheno_error))] <- 0
			if(is.null(ymin)){ymin <- min(c(pheno_means - pheno_error, 0), na.rm = TRUE)}
			if(is.null(ymax)){ymax <- max(pheno_means + pheno_error, na.rm = TRUE)}
			plot_height = ymax - ymin
				
			a <- barplot(pheno_means, ylim = c(ymin, ymax*1.1))
			segments(a, pheno_means-pheno_error, a, pheno_means+pheno_error, lwd = error_bar_lwd)
			abline(h = 0)
			#lower bar
			segments((a+error_bar_width), (pheno_means-pheno_error), (a-error_bar_width), pheno_means-pheno_error, lwd = error_bar_lwd)	
			#upper bar
			segments(a+error_bar_width, pheno_means+pheno_error, a-error_bar_width, pheno_means+pheno_error, lwd = error_bar_lwd)	
			par(xpd = TRUE)
			text(x = a[,1], y = ymin-(plot_height*0.05), labels = genotypes)
			mtext(pheno_name, side = 3, line = 1.5)

		} else { #instead if there are two markers

		
			#========================================================================
			# internal functions
			#========================================================================
	
			plot_bars_int <- function(pheno_vals, pheno_error, pheno_name, marker1, marker2, geno_n){
			
				if(length(pheno_vals) == 4){
					base_val = pheno_vals[1]
				}else{
					base_val = 0	
				}
					
				pred_add <- base_val + ((pheno_vals[2]-base_val) + (pheno_vals[3]-base_val))
				if(is.na(pred_add)){pred_add <- 0}
				pred_error <- pheno_error[2] + pheno_error[3]
		
				if(is.null(ymin)){
					if(error_bars != "none"){
						ymin <- min(c((pheno_vals-pheno_error), pred_add - pred_error, 0), na.rm = TRUE)
					}else{
						ymin <- min(c(pheno_vals, pred_add, 0), na.rm = TRUE)	
					}
				}
				if(is.null(ymax)){
					if(error_bars != "none"){
						ymax <- max(c((pheno_vals+pheno_error), pred_add + pred_error), na.rm = TRUE)
					}else{
						ymax <- max(c(pheno_vals, pred_add), na.rm = TRUE)	
					}
				}		
				
				
				full_range = ymax-ymin
				geno_n_y <- ymax*1.1
				ymax <- ymax*1.1
				
				label_y <- min(pheno_vals)-full_range*0.1
				a <- barplot(pheno_vals, ylim = c(ymin, ymax), axes = FALSE, 
				xlim = c(0, length(pheno_vals)*1.5), cex.names = 1.5, main = pheno_name, names.arg = NA)
				axis(2, cex.axis = 1.5)
				abline(h = 0)
				par(xpd = TRUE)
				text_cex = 1.5
				text(x = 0, y = ymin-(full_range*0.1), labels = marker1_label, adj = 1, cex = text_cex)
				text(x = 0, y = ymin-(full_range*0.22), labels = marker2_label, adj = 1, cex = text_cex)
				text(x = c(a), y = ymin-(full_range*0.1), labels = c(min_geno, max_geno, min_geno, max_geno), cex = text_cex)
				text(x = c(a), y = ymin-(full_range*0.22), labels = c(min_geno, min_geno, max_geno, max_geno), cex = text_cex)
				if(error_bars != "none"){
					segments(a, pheno_vals-pheno_error, a, pheno_vals+pheno_error, 
					lwd = error_bar_lwd)
					#lower bar
					segments((a+error_bar_width), (pheno_vals-pheno_error), 
					(a-error_bar_width), pheno_vals-pheno_error, lwd = error_bar_lwd)	
					#upper bar
					segments(a+error_bar_width, pheno_vals+pheno_error, a-error_bar_width,
					pheno_vals+pheno_error, lwd = error_bar_lwd)	
				}
				segments(a[length(a)]+addline_width, pred_add, a[length(a)]-addline_width, 
				pred_add, lty = 2, lwd = 2)
				poly_x <- c(a[length(a)]-addline_width, a[length(a)]+addline_width)
				poly_y <- c(pred_add-pred_error, pred_add+pred_error)
				polygon(x = c(poly_x, rev(poly_x)), y = rep(poly_y, each = 2), 
				col = rgb(253/256,192/256,134/256, alpha = 0.5))
				arrows(a[length(a)]+addline_width+addline_offset, pred_add, 
				a[length(a)]+addline_offset, pred_add, lty = 1, lwd = 2, length = 0.1)
				text(x = a[length(a)]+addline_width+addline_offset+text_offset, y = pred_add, 
				labels = "Additive\nPrediction", adj = 0)
				text(x = a[,1], y = rep(geno_n_y, length(a[,1])), labels = geno_n)
				par(xpd = FALSE)
				mtext(pheno_name, side = 2, line = 2.5)
			}			

		#========================================================================

		
		#the genotype combinations used are based on the
		#maximum genotype present. We want to show the 
		#maximum effect for the "mutant" phenotype.
		min_geno <- min(c(marker1_vals, marker2_vals), na.rm = TRUE)
		max_geno <- max(c(marker1_vals, marker2_vals), na.rm = TRUE)
				
				
		w_w <- intersect(which(marker1_vals == min_geno), which(marker2_vals == min_geno))
		w_m <- intersect(which(marker1_vals == max_geno), which(marker2_vals == min_geno))
		m_w <- intersect(which(marker1_vals == min_geno), which(marker2_vals == max_geno))
		m_m <- intersect(which(marker1_vals == max_geno), which(marker2_vals == max_geno))
				
		pheno_list <- vector(mode = "list", length = 4)
		names(pheno_list) <- c(paste(min_geno, min_geno, sep = "/"), 
				paste(min_geno, max_geno, sep = "/"), 
				paste(max_geno, min_geno, sep = "/"), 
				paste(max_geno, max_geno, sep = "/"))
	
		pheno_list[[1]] <- phenoV[w_w]
		pheno_list[[2]] <- phenoV[w_m]
		pheno_list[[3]] <- phenoV[m_w]
		pheno_list[[4]] <- phenoV[m_m]
											
		geno_n <- sapply(pheno_list, length)
		pheno_means <- sapply(pheno_list, function(x) mean(x, na.rm = TRUE))

		if(error_bars == "se"){
			pheno_error <- sapply(pheno_list, function(x) sd(x, na.rm = TRUE)/sqrt(length(x)))
		}
		if(error_bars == "sd"){
			pheno_error <- sapply(pheno_list, function(x) sd(x, na.rm = TRUE))
		}
		if(error_bars == "none"){
			pheno_error <- rep(0, length(pheno_list))
		}
		pheno_error[which(is.na(pheno_error))] <- 0

		if(!ref_centered){
			plot_bars_int(pheno_vals = pheno_means, pheno_error = pheno_error, 
			pheno_name = pheno_name, marker1 = marker1_vals, marker2 = marker2_vals, 
			geno_n = geno_n)
		}else{
			#also make a plot using B6 as the reference point
			#center on B6 genotype
			ref_pheno <- pheno_means-pheno_means[1]
			plot_bars_int(pheno_vals = ref_pheno, pheno_error, pheno_name = pheno_name, marker1 = marker1_vals, marker2 = marker2_vals, geno_n = geno_n)
		}		
	}
} #end case for two markers