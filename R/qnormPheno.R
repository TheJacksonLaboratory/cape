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
	
		
		layout.mat <- get.layout.mat(length(pheno.labels))
		layout(layout.mat)
		for(p in 1:length(pheno.which)){
			trait.locale <- which(colnames(data.obj$pheno) == pheno.which[p])
			trait.vals <- data.obj$pheno[,trait.locale]
			theo.norm <- rnorm(10000, mean(trait.vals, na.rm = TRUE), sd(trait.vals, na.rm = TRUE))
			qqplot(trait.vals, theo.norm, xlab = pheno.labels[p], ylab = "Theoretical Normal Quantiles",
			main = pheno.labels[p])
			abline(0,1)
			}
	}
