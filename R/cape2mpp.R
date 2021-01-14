#' Converts a \code{\link{read_population}} object to a multi-parent object
#'
#' This function converts an object formatted for cape 1.0
#' to an object formatted for cape 2.0
#'
#' @param data_obj a data_obj formatted for cape 1.0
#' @param geno_obj a genotype object. If geno_obj is NULL
#' the genotype object is generated from data_obj$geno.
#'
#' @return This function returns a list with two objects:
#' \code{list("data_obj" = data_obj, "geno_obj" = geno_obj)}
#' These two objects must be separated again to run through 
#' cape.
#'
#' @examples 
#' \dontrun{
#' new_data_obj <- cape2mpp(old_data_obj)
#' }
#' @export
cape2mpp <- function(data_obj, geno_obj = NULL){
	
	geno_locale <- which(names(data_obj) == "geno")
	if(length(geno_locale) > 0){
		geno <- data_obj$geno
		}else{
		geno <- geno_obj$geno	
		}
	
	geno_array <- array(NA, dim = c(nrow(geno),2,ncol(geno)))
	geno_array[,1,] <- 1-geno
	geno_array[,2,] <- geno	
	
	dimnames(geno_array) <- list(rownames(geno), c("A", "B"), colnames(geno))
	names(dimnames(geno_array)) <- c("mouse", "allele" ,"locus")
	
	geno_obj$geno <- geno_array
	
	data_obj$geno <- NULL
	data_obj$marker_names <- NULL

	geno_names <- list(rownames(data_obj$pheno), c("A", "B"), colnames(geno))
	names(geno_names) <- c("mouse", "allele", "locus")
	data_obj$geno_names <- geno_names


	results <- list("data_obj" = data_obj, "geno_obj" = geno_obj)
	return(results)

}