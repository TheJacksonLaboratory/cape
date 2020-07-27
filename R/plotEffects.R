#' Plot Interaction Effects
#'
#' This function plots phenotypic effects of 
#' individual cape interactions. It serves as
#' a wrapper for the functions \link{\code{plotLines}}
#' \link{\code{plotBars}} \link{\code{plotPoints}},
#' and \link{\code{IntHeat}}. Each of those functions
#' plots individual cape interactions in different forms.
#' 
#' @param data.obj A \link{\code{Cape}} object
#' @param geno.obj A genotype object
#' @param marker1 A string indicating the name 
#' of the source marker in the interaction. This can
#' also be the name of a covariate.
#' @param A string indicating the name of the target 
#' marker in the interaction. This can also be the name
#' of a covariate.
#' @param pheno.type One of "eigentraits", 
#' "normalized.traits", or "raw.traits", indicating which
#' traits to plot.
#' @param plot.type A letter referring to the desired style 
#' of the plot. The choices are the following: "l" - line plots,
#' "p" = points, "b" - bar plots, "h" - heat map.
#' @param error.bars The type of error bars to plot. Choices
#' are "none" (the default), "se" for standard error, or 
#' "sd" for standard deviation.
#' @param ymin A minimum value for the y axes across all plots.
#' If NULL, each y axis will be dertermined independently
#' @param ymax A maximum value for the y axes across all plots.
#' If NULL, each y axis will be dertermined independently
#' @param covar A vector of strings indicating which covariates,
#' if any, the traits should be adjusted for. If NULL, the 
#' covariates specified in the data.obj are used as default.
#' To prevent adjusting for covariates, use "none". 
#' @param marker1.label A string to use as the label for marker1
#' If NULL, the string used for marker1 will be used.
#' @param marker2.label A string to use as the label for marker2
#' If NULL, the string used for marker2 will be used.
#' @param bin.continuous.genotypes If TRUE, genotypes (and covariate)
#' values will be binned into 0, 0.5, and 1 values. This 
#' reduces the number of bins that traits need to be divided 
#' into, especially if there are only one or two individuals
#' with a 0.49 genotype, for example. Binning may not be
#' desirable when using the heatmap.
#' @param ref.centered A logical value indicating whether 
#' to center the values on the reference allele. Defaults 
#' to TRUE.
#' @param pheno A vector of strings to indicate which traits
#' to plot. If NULL, all traits are plotted.
#' @param gen.model One of "Additive", "Dominant", or "Recessive"
#' indicating how the genotypes should be coded. If Additive,
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
#' @param bins.marker1 Only used for heatmap plotting. The 
#' number of bins for marker1 if it is a continuously valued 
#' marker or covariate. The bins are used to fit a linear 
#' model and predict outcomes for a 2D grid of marker1 and
#' marker2 values. This argument can also be a vector of
#' bin values for binning at specific values.
#' @param bins.marker2 The same as bins.marker1, but for
#' marker2.
#' 
#' @details The "h" option calls \link{\code{IntHeat}}, which
#' fits linear models to each trait and both markers specified.
#' It uses those models to predict phenotype values along continuously
#' valued genotype bins and plots the predicted values as a heatmap.
#'
#' @return None
#'
#' @export

plotEffects <- function(data.obj, geno.obj, marker1, marker2 = NULL, 
pheno.type = "normalized", plot.type = c("l", "p", "b", "h"),
error.bars = "none", ymin = NULL, ymax = NULL, covar = NULL, 
marker1.label = NULL, marker2.label = NULL, bin.continuous.genotypes = TRUE, 
ref.centered = TRUE, gen.model = "Additive", bins.marker1 = 50, 
bins.marker2 = 50){

	plot.type = plot.type[1]
		
	#==========================================
	# get traits for plotting
	#==========================================
	covar.info <- get.covar(data.obj)		
	
	if(is.null(covar)){
		covar.names <- covar.info$covar.names
		}else{
		if(covar == "none"){
			covar.names = NULL
		}else{
			covar.names = covar
		}
	}
	
    pheno <- get.pheno(data.obj, pheno.type, covar.names)
	#==========================================


	#==========================================
	# get marker genotypes or covariate values 
	# for plotting
	#==========================================
	marker.vals <- get.marker.covar(data.obj, geno.obj, c(marker1, marker2))


   	#==========================================
	# line up genotypes and phenotypes
	#==========================================
    common.ind <- intersect(rownames(marker.vals), rownames(pheno))
    ind.pheno.locale <- match(common.ind, rownames(pheno))
    ind.geno.locale <- match(common.ind, rownames(marker.vals))
    
    pheno.to.plot <- pheno[ind.pheno.locale,,drop=FALSE]
    geno.to.plot <- marker.vals[ind.geno.locale,,drop=FALSE]
 	#==========================================


	#============================================================
	#bin the genotypes if necessary
	#============================================================
	if(bin.continuous.genotypes){
		geno.to.plot <- apply(geno.to.plot, 2, function(x) bin.vector(x, c(0, 0.5, 1)))
		}
	#============================================================


	#============================================================
	# Recode if specified
	#============================================================
		
	if(gen.model == "Dominant"){
		geno.to.plot[which(geno.to.plot >= 0.5)] <- 1
		geno.to.plot[which(geno.to.plot < 0.5)] <- 0
		}
	if(gen.model == "Recessive"){
		geno.to.plot[which(geno.to.plot <= 0.5)] <- 0			
		}
	#============================================================		


	#============================================================
	# assign the marker names
	#============================================================
	if(is.null(marker1.label)){marker1.label = marker1}
	if(is.null(marker2.label)){marker2.label = marker2}
	marker.names <- c(marker1.label, marker2.label)
	#============================================================

	
	#============================================================
	#figure out the layout for the plots
	#============================================================
	layout.mat <- get.layout.mat(ncol(pheno))
	layout(layout.mat)
	#============================================================
	

	for(ph in 1:ncol(pheno.to.plot)){
		phenoV = pheno.to.plot[,ph]
		pheno.name = colnames(pheno)[ph]
		marker1.vals <- geno.to.plot[,1]
		if(!is.null(marker2)){
			marker2.vals <- geno.to.plot[,2]
		}else{
			marker2.vals <- NULL
		}
		
		if(plot.type == "h"){
			if(is.null(marker2.vals)){stop("Two markers are required for the heat map.")}
			IntHeat(phenoV, marker1.vals, marker2.vals, pheno.name,
			marker1.label, marker2.label, bins1 = bins.marker1, bins2 = bins.marker2)
		}
		if(plot.type == "l"){
			plotLines(phenoV, marker1.vals, marker2.vals, pheno.name, 
			marker1.label, marker2.label, ymin, ymax, error.bars)
			}
		if(plot.type == "p"){
			plotPoints(phenoV, marker1.vals, marker2.vals, pheno.name, marker1.label,
			marker2.label, ymin, ymax)
			}
		if(plot.type == "b"){
			plotBars(phenoV, marker1.vals, marker2.vals, pheno.name, marker1.label,
			marker2.label, ymin, ymax, error.bars, ref.centered)
		}
		
	} #end looping through phenotypes


} #end function