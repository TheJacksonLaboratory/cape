#' converts a \code{\link{read.population}} object to a DO-enabled object
#'
#' This function converts an old cape object to
#' capeDO objects. The genotype object should ideally
#' be separate from the rest of the data
#'
#' @param data.obj a \code{\link{read.population}} object
#' @param geno.obj a genotype object. If this is not supplied then it is generated here.
#'
#' @return \code{list("data.obj" = data.obj, "geno.obj" = geno.obj)}
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