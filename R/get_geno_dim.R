#' Returns which dimensions the individual, locus, and
#' alleles are in in the genotype object.
#'
#' This is an internal function that returns the 
#' locations of genotype dimensions from the genotype
#' object. It is a relic from when DOQTL put loci in the
#' second dimension and alleles in the third, while R/qtl2
#' put loci in the third dimension and alleles in the second.
#'
#' @return Returns a vector of three named elements identifying
#' which dimension the individual (mouse), allele, and loci are in.
#' @keywords internal
#'

get_geno_dim <- function(){
	
	geno_dim <- c("mouse" = 1, "allele" = 2, "locus" = 3)
	return(geno_dim)
	
}