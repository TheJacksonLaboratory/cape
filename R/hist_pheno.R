#' Plot trait histograms
#' 
#' This function plots histograms of the traits held in 
#' the pheno slot of the data object.
#'
#' @param data_obj A \code{\link{Cape}} object
#' @param pheno_which A vector of strings indicating which 
#' traits to plot. Defaults to all traits.
#' @param pheno_labels A vector of strings providing alternate
#' names for the traits in the plot if the names in the data object
#' are not good for plotting
#'
#' @export

hist_pheno <- function(data_obj, pheno_which = NULL, pheno_labels = NULL){
	
	all_pheno <- data_obj$pheno

	if(is.null(pheno_which)){
		pheno_names <- colnames(all_pheno)
	}else{
		if(is.numeric(pheno_which)[1]){
			pheno_names <- colnames(all_pheno)[pheno_which]
		}else{
			pheno_names <- pheno_which	
		}
	}
	
	if(is.null(pheno_labels)){
		pheno_labels <- pheno_names
	}

	layout_mat <- get_layout_mat(length(pheno_names))
	layout(layout_mat)
	for(p in 1:length(pheno_labels)){
		hist(data_obj$pheno[,p], xlab = pheno_labels[p], main = pheno_labels[p])
	}
}
