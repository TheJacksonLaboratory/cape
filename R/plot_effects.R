#' Plot Interaction Effects
#'
#' This function plots phenotypic effects of 
#' individual cape interactions. It serves as
#' a wrapper for the functions \code{\link{plot_lines}}
#' \code{\link{plot_bars}} \code{\link{plot_points}},
#' and \code{\link{plot_int_heat}}. Each of those functions
#' plots individual cape interactions in different forms.
#' 
#' @param data_obj A \code{\link{Cape}} object
#' @param geno_obj A genotype object
#' @param marker1 A string indicating the name 
#' of the source marker in the interaction. This can
#' also be the name of a covariate.
#' @param marker2 Another string indicating the name 
#' of the source marker in the interaction. This can
#' also be the name of a covariate. Optional.
#' @param pheno_type One of "eigentraits", 
#' "normalized_traits", or "raw_traits", indicating which
#' traits to plot.
#' @param plot_type A letter referring to the desired style 
#' of the plot. The choices are the following: "l" - line plots,
#' "p" = points, "b" - bar plots, "h" - heat map.
#' @param error_bars The type of error bars to plot. Choices
#' are "none" (the default), "se" for standard error, or 
#' "sd" for standard deviation.
#' @param ymin A minimum value for the y axes across all plots.
#' If NULL, each y axis will be determined independently
#' @param ymax A maximum value for the y axes across all plots.
#' If NULL, each y axis will be dertermined independently
#' @param covar A vector of strings indicating which covariates,
#' if any, the traits should be adjusted for. If NULL, the 
#' covariates specified in the data_obj are used as default.
#' To prevent adjusting for covariates, use "none". 
#' @param marker1_label A string to use as the label for marker1
#' If NULL, the string used for marker1 will be used.
#' @param marker2_label A string to use as the label for marker2
#' If NULL, the string used for marker2 will be used.
#' @param bin_continuous_genotypes If TRUE, genotypes (and covariate)
#' values will be binned into 0, 0.5, and 1 values. This 
#' reduces the number of bins that traits need to be divided 
#' into, especially if there are only one or two individuals
#' with a 0.49 genotype, for example. Binning may not be
#' desirable when using the heatmap.
#' @param ref_centered A logical value indicating whether 
#' to center the values on the reference allele. Defaults 
#' to TRUE.
#' @param gen_model1 One of "Additive", "Dominant", or "Recessive"
#' indicating how the genotype should be coded for the first 
#' marker. If Additive,
#' genotypes are coded as 0 for homozygous reference allele,
#' 1 for homozygous alternate allele, and 0.5 for heterozygous.
#' If Dominant, any allele probability greater than 0.5 is 
#' set to 1. If recessive, any allele probability less than
#' or equal to 0.5 is set to 0. In other words, for the 
#' dominant coding, heterozygotes are grouped with the 
#' homozygous alternate genotypes: 0 vs. (0.5,1). This shows
#' the effect of having any dose of the alternate allele. With
#' a recessive coding, heterozygotes are grouped with the
#' homozygous reference genotypes: (0, 0.5) vs. 1. This shows
#' the effect of having two copies of the alternate allele 
#' vs. having fewer than two copies.
#' @param gen_model2 The same as gen_model1, but for the second 
#' marker.
#' @param bins_marker1 Only used for heatmap plotting. The 
#' number of bins for marker1 if it is a continuously valued 
#' marker or covariate. The bins are used to fit a linear 
#' model and predict outcomes for a 2D grid of marker1 and
#' marker2 values. This argument can also be a vector of
#' bin values for binning at specific values.
#' @param bins_marker2 The same as bins_marker1, but for
#' marker2.
#' 
#' @details The "h" option calls \code{\link{plot_int_heat}}, which
#' fits linear models to each trait and both markers specified.
#' It uses those models to predict phenotype values along continuously
#' valued genotype bins and plots the predicted values as a heatmap.
#'
#' @return None
#' 
#' @examples 
#' \dontrun{
#' marker1 <- dimnames(geno_obj)[[3]][1]
#' marker2 <- dimnames(geno_obj)[[3]][2]
#' plot_effects(data_obj, geno_obj, plot_type = "l", error_bars = "se")
#' }
#'
#' @export

plot_effects <- function(data_obj, geno_obj, marker1, marker2 = NULL, 
	pheno_type = "normalized", plot_type = c("l", "p", "b", "h"),
	error_bars = "none", ymin = NULL, ymax = NULL, covar = NULL, 
	marker1_label = NULL, marker2_label = NULL, bin_continuous_genotypes = TRUE, 
	ref_centered = TRUE, gen_model1 = "Additive", gen_model2 = "Additive",
	bins_marker1 = 50, bins_marker2 = 50){

	plot_type = plot_type[1]
		
	#==========================================
	# get traits for plotting
	#==========================================
	covar_info <- get_covar(data_obj)		
	
	if(is.null(covar)){
		covar_names <- covar_info$covar_names
	}else{
		if(covar == "none"){
			covar_names = NULL
		}else{
			covar_names = covar
		}
	}
	
    pheno <- get_pheno(data_obj, pheno_type, covar_names)
	#==========================================


	#==========================================
	# get marker genotypes or covariate values 
	# for plotting
	#==========================================
	marker_vals <- get_marker_covar(data_obj, geno_obj, c(marker1, marker2))


   	#==========================================
	# line up genotypes and phenotypes
	#==========================================
    common_ind <- intersect(rownames(marker_vals), rownames(pheno))
    ind_pheno_locale <- match(common_ind, rownames(pheno))
    ind_geno_locale <- match(common_ind, rownames(marker_vals))
    
    pheno_to_plot <- pheno[ind_pheno_locale,,drop=FALSE]
    geno_to_plot <- marker_vals[ind_geno_locale,,drop=FALSE]
 	#==========================================


	#============================================================
	#bin the genotypes if necessary
	#============================================================
	if(bin_continuous_genotypes){
		geno_to_plot <- apply(geno_to_plot, 2, function(x) bin_vector(x, c(0, 0.5, 1)))
	}
	#============================================================


	#============================================================
	# Recode if specified
	#============================================================
		
	if(gen_model1 == "Dominant"){
		geno_to_plot[which(geno_to_plot[,1] >= 0.5),1] <- 1
		geno_to_plot[which(geno_to_plot[,1] < 0.5),1] <- 0
	}
	if(gen_model1 == "Recessive"){
		geno_to_plot[which(geno_to_plot[,1] <= 0.5),1] <- 0			
	}

	if(gen_model2 == "Dominant"){
		geno_to_plot[which(geno_to_plot[,2] >= 0.5),2] <- 1
		geno_to_plot[which(geno_to_plot[,2] < 0.5),2] <- 0
	}
	if(gen_model2 == "Recessive"){
		geno_to_plot[which(geno_to_plot[,2] <= 0.5),2] <- 0			
	}

	#============================================================		



	#============================================================
	# assign the marker names
	#============================================================
	if(is.null(marker1_label)){marker1_label = marker1}
	if(is.null(marker2_label)){marker2_label = marker2}
	marker_names <- c(marker1_label, marker2_label)
	#============================================================

	
	#============================================================
	#figure out the layout for the plots
	#============================================================
	layout_mat <- get_layout_mat(ncol(pheno))
	layout(layout_mat)
	#============================================================
	

	for(ph in 1:ncol(pheno_to_plot)){
		phenoV = pheno_to_plot[,ph]
		pheno_name = colnames(pheno)[ph]
		marker1_vals <- geno_to_plot[,1]
		if(!is.null(marker2)){
			marker2_vals <- geno_to_plot[,2]
		}else{
			marker2_vals <- NULL
		}
		
		if(plot_type == "h"){
			if(is.null(marker2_vals)){stop("Two markers are required for the heat map.")}
		  plot_int_heat(phenoV, marker1_vals, marker2_vals, pheno_name,
			marker1_label, marker2_label, bins1 = bins_marker1, bins2 = bins_marker2)
		}
		if(plot_type == "l"){
		  plot_lines(phenoV, marker1_vals, marker2_vals, pheno_name, 
			marker1_label, marker2_label, ymin, ymax, error_bars)
		}
		if(plot_type == "p"){
		  plot_points(phenoV, marker1_vals, marker2_vals, pheno_name, marker1_label,
			marker2_label, ymin, ymax)
		}
		if(plot_type == "b"){
		  plot_bars(phenoV, marker1_vals, marker2_vals, pheno_name, marker1_label,
			marker2_label, ymin, ymax, error_bars, ref_centered)
		}
		
	} #end looping through phenotypes


} #end function