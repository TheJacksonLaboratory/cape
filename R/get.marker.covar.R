#' Get genotype or covariate values
#' 
#' In some cape results plotting functions
#' either marker genotypes or covariates can
#' be used for plotting effects on traits. 
#' However markers and covariates are stored in
#' different places so it can be annoying to 
#' retrieve their values.
#' 
#' This function looks both in the genotype
#' matrix and in covariate tables to find those
#' values.
#' 
#' @param data.obj A \code{\link{Cape}} object
#' @param geno.obj A genotype object
#' @param marker.covar.names A vector of strings 
#' which can contain marker names with alleles appended
#' or covariate names. 
#' 
#' @return This function returns a matrix with individuals
#' in rows and markers/covariates in columns.
#'

 
get.marker.covar <- function(data.obj, geno.obj, marker.covar.names){

	geno.names <- dimnames(geno.obj)

	split.markers <- strsplit(marker.covar.names, "_")
	just.markers <- sapply(split.markers, function(x) x[1])
	just.alleles <- sapply(split.markers, function(x) x[2])
    marker.locale <- match(just.markers, geno.names[[3]])
    allele.locale <- match(just.alleles, geno.names[[2]])

	#================================================
	# align individuals from genotype and covariate
	# matrices
	#================================================
	covar.info <- get.covar(data.obj)
	common.ind <- intersect(rownames(covar.info$covar.table), rownames(geno.obj))
	covar.ind.locale <- match(common.ind, rownames(covar.info$covar.table))
	geno.ind.locale <- match(common.ind, rownames(geno.obj))

	#================================================
	# get genotype values for all genetic markers
	#================================================
    marker.vals <- sapply(1:length(just.markers), function(x) geno.obj[geno.ind.locale,allele.locale[x], marker.locale[x]])


	#================================================
	#check for any markers without alleles specified
	#================================================
	na.markers <- which(is.na(marker.locale))
	na.alleles <- which(is.na(allele.locale))
	marker.no.allele <- setdiff(na.alleles, na.markers)
	#if there are markers without alleles set,
	#set allele locale to 2 with a message
	if(length(marker.no.allele) > 0){ 
		cat('Setting missing alleles to', geno.names[[2]][2], "\n")
		allele.locale[marker.no.allele] <- 2
	}

    
	#================================================
	# for any markers whose position couldn't be 
	# found look among the covariates.
	#================================================
	
    for(i in na.markers){
	  	covar.locale <- which(covar.info$covar.names == just.markers[i])
	  	if(is.na(covar.locale)){
		  	cat("Cannot find", marker.covar.names[i], "\n")
		  	}else{
		  	marker.vals[,i] <- covar.info$covar.table[covar.ind.locale,covar.locale]
		  	}
	  	}
    
	colnames(marker.vals) <- marker.covar.names
	rownames(marker.vals) <- common.ind
	return(marker.vals)
}