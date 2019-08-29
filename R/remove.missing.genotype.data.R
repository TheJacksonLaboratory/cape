
#data.obj = cross; geno.obj = geno; ind.missing.thresh = 0; marker.missing.thresh = 0; prioritize = c("ind", "marker", "fewer")

#' Removes missing marker (column) or individual (row) data
#' 
#' This function removes missing data and let you remove 
#' individuals before markers, vice versa, or 
#' decide based on which will remove fewer elements.
#' because get.geno uses geno.names in data.obj to
#' retrieve the genotype.object, nothing needs to be
#' removed from the geno.obj
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#' @param ind.missing.thresh percent of individuals that are acceptible to remove, default = 0
#' @param marker.missing.thresh percent of genotype markers that are acceptible to remove, default = 0
#' @param prioritize the basis prioritization is one of 
#'        "fewer" = remove the fewest possible cells from the matrix
#'        "ind" = remove the fewest possible individuals
#'        "marker" = remove the fewest possible markers
#'
#' @return an updated cape object
#'
#' @export
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