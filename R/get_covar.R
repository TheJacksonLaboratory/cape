#' Get covariate information
#' 
#' This function returns information about the covariates
#' specified for the cape run.
#'
#' @param data_obj a \code{\link{Cape}} object
#'
#' @return Returns a list with the following elements:
#' covar_names: a character vector holding the names of the covariates  
#' covar_type: a character vector indicating whether each covariate
#' derived from the phenotype matrix ("p") or the genotype matrix ("g")
#' covar_loc: A numeric vector indicating the locations of each covariate
#' covar_table: A matrix holding the individual values for each covariate.
#' 
#' @export

get_covar <- function(data_obj){
		
	covar_table <- cbind(data_obj$p_covar_table, data_obj$g_covar_table)
	
	if(!is.null(data_obj$p_covar_table)){
		num_p_covar <- dim(data_obj$p_covar_table)
		if (length(num_p_covar) == 1) {
		  # there's only one covar column
		  num_p_covar <- 1
		} else {
		  num_p_covar <- num_p_covar[2]
		}
		}else{
		num_p_covar <- 0	
		}
		
	if(!is.null(data_obj$g_covar_table)){
		num_g_covar <- dim(data_obj$g_covar_table)[2]
		}else{
		num_g_covar <- 0	
		}
	
	covar_type <- c(rep("p", num_p_covar), rep("g", num_g_covar))

	covar <- as.vector(c(data_obj$p_covar, data_obj$g_covar[1,]))
	
		covar_names <- c(data_obj$p_covar, data_obj$g_covar[1,])
		p_covar_loc <- NULL
		g_covar_loc <- NULL
				
		p_covar_locale <- which(covar_type == "p")
		g_covar_locale <- which(covar_type == "g")
		
		if(length(p_covar_locale) > 0){
			p_covar_loc <- rep(1:length(which(covar_type == "p")))
			}
		if(length(g_covar_locale) > 0){
			g_covar_loc <- data_obj$g_covar["position",]
			}
		covar_loc <- c(p_covar_loc, g_covar_loc)
		covar_table <- cbind(data_obj$p_covar_table, data_obj$g_covar_table)
		
		if(length(covar_loc) < length(covar)){
			not_found <- setdiff(covar, covar_names)
			if(length(not_found) > 0){
				warning("I could not find the following covariates: ", paste(not_found, collapse = ", "))
				stop()
				}
			}
	
	result <- list("covar_names" = covar_names, "covar_type" = covar_type, "covar_loc" = covar_loc, "covar_table" = covar_table)
	return(result)	
	
}