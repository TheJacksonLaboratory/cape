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
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#' @param ind.missing.thresh Allowable amount of missing information for an individual. 
#' If 10%, only individuals missing more than 10% of markers will be removed. If 0%, the
#' default, all individuals with any missing data at all will be removed.
#' @param marker.missing.thresh Allowable amount of missing information for a marker. 
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
#' the function \link{\code{get.geno}} should return an array with no missing data if ind.missing.thresh
#' and marker.missing.thresh are both 0. If these numbers are higher, no individual or marker will
#' be missing more than the set percentage of data.
#' 
#' details All missing genotype data must either be imputed or removed if using the kinship correction.
#' Running \link{\code{impute.missing.geno}} prior to running \link{\code{remove.missing.genotype.data}}
#' ensures that the least possible amount of data are removed before running cape. In some cases, there
#' will be missing genotype data even after running \link{\code{impute.missing.geno}}, in which case,
#' \code{remove.missing.genotype.data} still needs to be run. 
#' The function \link{\code{run.cape}} automatically runs these steps when \code{use_kinship}
#' is set to TRUE. 
#' 
#' @seealso \link{\code{get.geno}}, \link{\code{impute.missing.geno}}, \link{\code{run.cape}}
#' @export
#'
remove.missing.genotype.data <- function(data.obj, geno.obj = NULL, ind.missing.thresh = 0,
                                         marker.missing.thresh = 0, prioritize = c("fewer", "ind", "marker")){
	
	geno <- get.geno(data.obj, geno.obj)
	num.na <- length(which(is.na(geno)))
	if(num.na == 0){cat("No missing genotypes\n"); return(data.obj)}
		
	prioritize <- prioritize[1]
	
	#========================================================================
	# internal functions
	#========================================================================
	assess.missing <- function(data.obj, geno.obj){
		geno <- get.geno(data.obj, geno.obj)

		flat.geno <- flatten.array(arrayX = geno, 3, 1, "mean")
		
		ind.missing.percent <- as.vector(apply(flat.geno, 1, function(x) length(which(is.na(x)))))/dim(flat.geno)[2]*100
		marker.missing.percent <- as.vector(apply(flat.geno, 2, function(x) length(which(is.na(x)))))/dim(flat.geno)[1]*100
		
		ind.missing.lots <- which(ind.missing.percent > ind.missing.thresh)
		marker.missing.lots <- which(marker.missing.percent > marker.missing.thresh)
		
		results <- list(ind.missing.lots, marker.missing.lots)
		names(results) <- c("ind.missing.lots", "marker.missing.lots")
		return(results)
		}	
	
		#========================================================================
		# end internal functions
		#========================================================================
		
	
		#========================================================================
		# If we are removing individuals before markers...
		#========================================================================
		if(length(grep("i", prioritize, ignore.case = TRUE)) > 0){
			test <- assess.missing(data.obj, geno)
			if(length(test$ind.missing.lots) > 0){
			cat(paste("Removing ", length(test$ind.missing.lots), " individual(s) with more than ", ind.missing.thresh, "% missing data.\n", sep = ""))
				data.obj <- remove.ind(data.obj, ind.to.remove = test$ind.missing.lots)
				}
			test <- assess.missing(data.obj, geno)
			if(length(test$marker.missing.lots) > 0){
				cat(paste("Removing ", length(test$marker.missing.lots), " markers with more than ", marker.missing.thresh, "% missing data.\n", sep = ""))		
				data.obj <- remove.markers(data.obj, markers.to.remove = test$marker.missing.lots)
				}
			}


		#========================================================================
		# If we are removing markers before individuals...
		#========================================================================	
		if(length(grep("m", prioritize, ignore.case = TRUE)) > 0){
			test <- assess.missing(data.obj, geno)
			if(length(test$marker.missing.lots) > 0){
				cat(paste("Removing ", length(test$marker.missing.lots), " markers with more than ", marker.missing.thresh, "% missing data.\n", sep = ""))		
				data.obj <- remove.markers(data.obj, markers.to.remove = test$marker.missing.lots)
				}
			test <- assess.missing(data.obj, geno)
			if(length(test$ind.missing.lots) > 0){
			cat(paste("Removing ", length(test$ind.missing.lots), " individual(s) with more than ", ind.missing.thresh, "% missing data.\n", sep = ""))
				data.obj <- remove.ind(data.obj, ind.to.remove = test$ind.missing.lots)
				}
			}
			
			
		#========================================================================
		# If we are removing whichever has fewer missing
		#========================================================================			
		if(length(grep("f", prioritize, ignore.case = TRUE)) > 0){
			test <- assess.missing(data.obj, geno)
			perc.ind <- length(test$ind.missing.lots)/nrow(cross$pheno)
			perc.markers <- length(test$marker.missing.lots)/length(data.obj$geno_names[[3]])
			
			if(perc.markers > 0 || perc.ind > 0){
				if(perc.markers < perc.ind){
					cat(paste("Removing ", length(test$marker.missing.lots), " markers with more than ", marker.missing.thresh, "% missing data.\n", sep = ""))		
					data.obj <- remove.markers(data.obj, markers.to.remove = test$marker.missing.lots)
					}else{
				cat(paste("Removing ", length(test$ind.missing.lots), " individual(s) with more than ", ind.missing.thresh, "% missing data.\n", sep = ""))
					data.obj <- remove.ind(data.obj, ind.to.remove = test$ind.missing.lots)
					}
				}
			}

			return(data.obj)

}