#' Plot trait histograms
#' 
#' This function plots histograms of the traits held in 
#' the pheno slot of the data object.
#'
#' @param data.obj A \code{\link{Cape}} object
#' @param pheno.which A vector of strings indicating which 
#' traits to plot. Defaults to all traits.
#' @param pheno.labels A vector of strings providing alternate
#' names for the traits in the plot if the names in the data object
#' are not good for plotting
#'
#' @export

histPheno <- function(data.obj, pheno.which = NULL, pheno.labels = NULL){
	
	all.pheno <- data.obj$pheno


	if(is.null(pheno.which)){
		pheno.names <- colnames(all.pheno)
		}else{
			if(is.numeric(pheno.which)[1]){
				pheno.names <- colnames(all.pheno)[pheno.which]
				}else{
				pheno.names <- pheno.which	
				}
			}
	
	if(is.null(pheno.labels)){
		pheno.labels <- pheno.names
		}
		

		layout.mat <- get.layout.mat(length(pheno.names))
		layout(layout.mat)
		for(p in 1:length(pheno.labels)){
			hist(data.obj$pheno[,p], xlab = pheno.labels[p], main = pheno.labels[p])
			}
	}
