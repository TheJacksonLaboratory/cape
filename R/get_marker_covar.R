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
#' @param data_obj A \code{\link{Cape}} object
#' @param geno_obj A genotype object
#' @param marker_covar_names A vector of strings 
#' which can contain marker names with alleles appended
#' or covariate names. 
#' 
#' @return This function returns a matrix with individuals
#' in rows and markers/covariates in columns.
#' @keywords internal
#'

 
get_marker_covar <- function(data_obj, geno_obj, marker_covar_names){

	geno_names <- dimnames(geno_obj)

	split_markers <- strsplit(marker_covar_names, "_")
	just_markers <- sapply(split_markers, function(x) x[1])
	just_alleles <- sapply(split_markers, function(x) x[2])
    marker_locale <- match(just_markers, geno_names[[3]])
    allele_locale <- match(just_alleles, geno_names[[2]])

	#================================================
	# align individuals from genotype and phenotype
	# matrices
	#================================================
	common_ind <- intersect(rownames(data_obj$pheno), rownames(geno_obj))
	geno_ind_locale <- match(common_ind, rownames(geno_obj))

	marker_vals <- sapply(1:length(just_markers), function(x) geno_obj[geno_ind_locale,allele_locale[x], marker_locale[x]])

	#================================================
	# get genotype values for all genetic markers
	#================================================
	if(all(!is.na(marker_locale))){
		return(marker_vals)
	}else{ 

		#================================================
		# if some values are NAs, these are probably 
		# covariates
		#================================================

		na_markers <- which(is.na(marker_locale))
	    
		#================================================
		# for any markers whose position couldn't be 
		# found look among the covariates.
		#================================================
		covar_info <- get_covar(data_obj)
		covar_ind_locale <- match(common_ind, rownames(covar_info$covar_table))
		for(i in na_markers){
				covar_locale <- which(covar_info$covar_names == just_markers[i])
				if(is.na(covar_locale)){
					warning("Cannot find", marker_covar_names[i], "\n")
				}else{
					marker_vals[,i] <- covar_info$covar_table[covar_ind_locale,covar_locale]
				}
			}
		colnames(marker_vals) <- marker_covar_names
		rownames(marker_vals) <- common_ind
		return(marker_vals)
	}
}