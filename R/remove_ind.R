#' Remove individuals
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param ind_to_remove Indices of individuals to remove
#' @param names_to_remove Names of individuals to remove
#' Only one of ind_to_remove or names_to_remove should be specified.
#'
#' @return an updated cape data object with specified individuals removed.
#' 
#' @examples 
#' \dontrun{
#' #remove males
#' covar_info <- get_covar(data_obj)
#' male_idx <- which(covar_info$covar_table[,"sex"] == 1)
#' data_obj <- remove_ind(data_obj, ind_to_remove = male_idx)
#' }
#'
#' @export
remove_ind <- function(data_obj, ind_to_remove = NULL, names_to_remove = NULL){

	ind_idx <- NULL

	if(!is.null(names_to_remove)){
		ind_idx <- match(names_to_remove, rownames(data_obj$pheno))
	}
	if(!is.null(ind_to_remove)){
		ind_idx = unique(c(ind_idx, ind_to_remove))
	}

	if(length(ind_idx) > 0){
		#remove individuals from the phenotype matrix in data_obj
		#and geno_names
		data_obj$pheno <- data_obj$pheno[-ind_idx,,drop=FALSE]
		data_obj$geno_names[[1]] <- data_obj$geno_names[[1]][-ind_idx]

		#if covariates have already been assigned, remove individuals
		#from these tables as well.
		if(!is.null(data_obj$p_covar_table)){
			tryCatch(
			{
				data_obj$p_covar_table <- data_obj$p_covar_table[-ind_idx,,drop=FALSE]
			},
			error=function(cond) {
				# if the table is down to one column we have to slice it differently
				data_obj$p_covar_table <- data_obj$p_covar_table[-ind_idx,drop=FALSE]
			}
			)
		}
		if(!is.null(data_obj$g_covar_table)){
			data_obj$g_covar_table <- data_obj$g_covar_table[-ind_idx,,drop=FALSE]
		}
		if(!is.null(data_obj$raw_pheno)){
			data_obj$raw_pheno <- data_obj$raw_pheno[-ind_idx,,drop=FALSE]
		}
		if(!is.null(data_obj$ET)){
			warning("get_eigentraits needs to be re-run because individuals were removed.\n")
		}
	}
	return(data_obj)
}