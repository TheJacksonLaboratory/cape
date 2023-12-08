#' Plot trait distributions
#' 
#' This function plots the quantiles of each trait
#' against quantiles of a theoretical normal distribution.
#' This provides a way to check whether traits are normally
#' distributed
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param pheno_which A vector of trait names to plot. The default is to plot all traits.
#' @param pheno_labels A vector of names for traits to appear in the plot in case the column names are not very pretty.
#' 
#' @importFrom graphics abline
#' @importFrom stats rnorm qqplot
#' 
#' @export

qnorm_pheno <- function(data_obj, pheno_which = NULL, pheno_labels = NULL){
	
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
	
		
	layout_mat <- get_layout_mat(length(pheno_labels))
	layout(layout_mat)
	for(p in 1:length(pheno_which)){
		trait_locale <- which(colnames(data_obj$pheno) == pheno_which[p])
		trait_vals <- data_obj$pheno[,trait_locale]
		theo_norm <- rnorm(10000, mean(trait_vals, na.rm = TRUE), sd(trait_vals, na.rm = TRUE))
		qqplot(trait_vals, theo_norm, xlab = pheno_labels[p], ylab = "Theoretical Normal Quantiles",
		main = pheno_labels[p])
		abline(0,1)
	}
}
