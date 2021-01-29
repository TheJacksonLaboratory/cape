#' Removes individuals and/or markers with missing data
#' 
#' Because there an be no missing data when calculating
#' the kinship correction, we need a way to remove either
#' individuals or markers with missing data. We also need
#' a way to calculate which of these options will remove the
#' least amount of data.
#'
#' For example, if there is one marker with no data at all,
#' we would rather remove that one marker, than all individuals
#' with missing data. Alternatively, if there is one individual
#' with very sparse genotyping, we would prefer to remove that 
#' single individual, rather than all markers with missing data.
#' 
#' This function provides a way to calculate whether individuals
#' or markers should be prioritized when removing data. It then
#' removes those individuals or markers.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' @param ind_missing_thresh Allowable amount of missing information for an individual. 
#' If 10%, only individuals missing more than 10% of markers will be removed. If 0%, the
#' default, all individuals with any missing data at all will be removed.
#' @param marker_missing_thresh Allowable amount of missing information for a marker. 
#' If 10%, only markers missing more than 10% of of genotypes will be removed. If 0%, the
#' default, all markers with any missing data at all will be removed.
#' @param prioritize the basis prioritization is one of 
#'        "fewer" = calculate whether removing individuals or markers 
#'           will remove fewer data points, and start with that.
#'        "ind" = remove individuals with missing data before considering 
#'           markers with missing data.
#'        "marker" = remove markers with missing data before considering individuals.
#'
#' @return The cape object is returned with individuals and markers removed. After this step,
#' the function \code{\link{get_geno}} should return an array with no missing data if ind_missing_thresh
#' and marker_missing_thresh are both 0. If these numbers are higher, no individual or marker will
#' be missing more than the set percentage of data.
#' 
#' details All missing genotype data must either be imputed or removed if using the kinship correction.
#' Running \code{\link{impute_missing_geno}} prior to running \code{\link{remove_missing_genotype_data}}
#' ensures that the least possible amount of data are removed before running cape. In some cases, there
#' will be missing genotype data even after running \code{\link{impute_missing_geno}}, in which case,
#' \code{remove_missing_genotype_data} still needs to be run. 
#' The function \code{\link{run_cape}} automatically runs these steps when \code{use_kinship}
#' is set to TRUE. 
#' 
#' @seealso \code{\link{get_geno}}, \code{\link{impute_missing_geno}}, \code{\link{run_cape}}
#' 
#' @examples 
#' \dontrun{
#' #remove entries with more than 10\% missing data prioitizing 
#' #removal of markers
#' data_obj <- remove_missing_genotype_data(data_obj, geno_obj, 
#' marker_missing_thresh = 10, ind_missing_thresh = 10,
#' prioritize = "marker")
#' 
#' #remove markers with more than 5\% missing data and markers with 
#' #more than 50\% 
#' #missing data, prioritizing removal of individuals.
#' data_obj <- remove_missing_genotype_data(data_obj, geno_obj, 
#' ind_missing_thresh = 10, marker_missing_thresh = 50,
#' prioritize = "ind")
#' 
#' #remove entries witn any missing data prioritizing whichever 
#' #method removes the least amount of data
#' data_obj <- remove_missing_genotype_data(data_obj, geno_obj)
#' 
#' }
#' 
#' @export
remove_missing_genotype_data <- function(data_obj, geno_obj = NULL, ind_missing_thresh = 0,
                                         marker_missing_thresh = 0, prioritize = c("fewer", "ind", "marker")){
	
	geno <- get_geno(data_obj, geno_obj)
	num_na <- length(which(is.na(geno)))
	if(num_na == 0){message("No missing genotypes\n"); return(data_obj)}
		
	prioritize <- prioritize[1]
	
	#========================================================================
	# internal functions
	#========================================================================
	assess_missing <- function(data_obj, geno_obj){
		geno <- get_geno(data_obj, geno_obj)

		flat_geno <- flatten_array(arrayX = geno, 3, 1, "mean")
		
		ind_missing_percent <- as.vector(apply(flat_geno, 1, function(x) length(which(is.na(x)))))/dim(flat_geno)[2]*100
		marker_missing_percent <- as.vector(apply(flat_geno, 2, function(x) length(which(is.na(x)))))/dim(flat_geno)[1]*100
		
		ind_missing_lots <- which(ind_missing_percent > ind_missing_thresh)
		marker_missing_lots <- which(marker_missing_percent > marker_missing_thresh)
		
		results <- list(ind_missing_lots, marker_missing_lots)
		names(results) <- c("ind_missing_lots", "marker_missing_lots")
		return(results)
	}
	
	#========================================================================
	# end internal functions
	#========================================================================
		
	
	#========================================================================
	# If we are removing individuals before markers...
	#========================================================================
	if(length(grep("i", prioritize, ignore.case = TRUE)) > 0){
		test <- assess_missing(data_obj, geno)
		if(length(test$ind_missing_lots) > 0){
			message(paste("Removing ", length(test$ind_missing_lots), " individual(s) with more than ", ind_missing_thresh, "% missing data.\n", sep = ""))
			data_obj <- remove_ind(data_obj, ind_to_remove = test$ind_missing_lots)
		}
		test <- assess_missing(data_obj, geno)
		if(length(test$marker_missing_lots) > 0){
			message(paste("Removing ", length(test$marker_missing_lots), " markers with more than ", marker_missing_thresh, "% missing data.\n", sep = ""))
			data_obj <- remove_markers(data_obj, markers_to_remove = test$marker_missing_lots)
		}
	}


	#========================================================================
	# If we are removing markers before individuals...
	#========================================================================	
	if(length(grep("m", prioritize, ignore.case = TRUE)) > 0){
		test <- assess_missing(data_obj, geno)
		if(length(test$marker_missing_lots) > 0){
			message(paste("Removing ", length(test$marker_missing_lots), " markers with more than ", marker_missing_thresh, "% missing data.\n", sep = ""))		
			data_obj <- remove_markers(data_obj, markers_to_remove = test$marker_missing_lots)
		}
		test <- assess_missing(data_obj, geno)
		if(length(test$ind_missing_lots) > 0){
			message(paste("Removing ", length(test$ind_missing_lots), " individual(s) with more than ", ind_missing_thresh, "% missing data.\n", sep = ""))
			data_obj <- remove_ind(data_obj, ind_to_remove = test$ind_missing_lots)
		}
	}
		
		
	#========================================================================
	# If we are removing whichever has fewer missing
	#=======================================================================
	if(length(grep("f", prioritize, ignore.case = TRUE)) > 0){
		test <- assess_missing(data_obj, geno)
		perc_ind <- length(test$ind_missing_lots)/nrow(data_obj$pheno)
		perc_markers <- length(test$marker_missing_lots)/length(data_obj$geno_names[[3]])
		
		if(perc_markers > 0 || perc_ind > 0){
			if(perc_markers < perc_ind){
				message(paste("Removing ", length(test$marker_missing_lots), " markers with more than ", marker_missing_thresh, "% missing data.\n", sep = ""))		
				data_obj <- remove_markers(data_obj, markers_to_remove = test$marker_missing_lots)
			}else{
				message(paste("Removing ", length(test$ind_missing_lots), " individual(s) with more than ", ind_missing_thresh, "% missing data.\n", sep = ""))
				data_obj <- remove_ind(data_obj, ind_to_remove = test$ind_missing_lots)
			}
		}
	}

	return(data_obj)

}