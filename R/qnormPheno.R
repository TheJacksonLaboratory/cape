#' Plot trait distributions
#' 
#' This function plots the quantiles of each trait
#' against quantiles of a theoretical normal distribution.
#' This provides a way to check whether traits are normally
#' distributed
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param pheno.which A vector of trait names to plot. The default is to plot all traits.
#' @param pheno.labels A vector of names for traits to appear in the plot in case the column names are not very pretty.
#' 
#' @export

qnormPheno <- function(data.obj, pheno.which = NULL, pheno.labels = NULL){
	
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
		
	pheno.num.pairs <- pair.matrix(1:length(pheno.labels))

		layout.mat <- get.layout.mat(dim(pheno.num.pairs)[1])
		layout(layout.mat)
		for(p in 1:nrow(pheno.num.pairs)){
			qqplot(data.obj$pheno[,pheno.num.pairs[p,1]], data.obj$pheno[,pheno.num.pairs[p,2]], xlab = pheno.labels[pheno.num.pairs[p,1]], ylab = pheno.labels[pheno.num.pairs[p,2]])
			}
	}
