#' Converts a \code{\link{read.population}} object to a multi-parent object
#'
#' This function converts an object formatted for cape 1.0
#' to an object formatted for cape 2.0
#'
#' @param data.obj a data.obj formatted for cape 1.0
#' @param geno.obj a genotype object. If geno.obj is NULL
#' the genotype object is generated from data.obj$geno.
#'
#' @return This function returns a list with two objects:
#' \code{list("data.obj" = data.obj, "geno.obj" = geno.obj)}
#' These two objects must be separated again to run through 
#' cape.
#'
#' @export
cape2mpp <- function(data.obj, geno.obj = NULL){
	
	geno.locale <- which(names(data.obj) == "geno")
	if(length(geno.locale) > 0){
		geno <- data.obj$geno
		}else{
		geno <- geno.obj$geno	
		}
	
	geno.array <- array(NA, dim = c(nrow(geno),2,ncol(geno)))
	geno.array[,1,] <- 1-geno
	geno.array[,2,] <- geno	
	
	dimnames(geno.array) <- list(rownames(geno), c("A", "B"), colnames(geno))
	names(dimnames(geno.array)) <- c("mouse", "allele" ,"locus")
	
	geno.obj$geno <- geno.array
	
	data.obj$geno <- NULL
	data.obj$marker_names <- NULL

	geno_names <- list(rownames(data.obj$pheno), c("A", "B"), colnames(geno))
	names(geno_names) <- c("mouse", "allele", "locus")
	data.obj$geno_names <- geno_names


	results <- list("data.obj" = data.obj, "geno.obj" = geno.obj)
	return(results)

}