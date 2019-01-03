#This function writes the parameters 
#used in a cape analysis

readParameters <- function(filename = "cape.parameters.txt"){

	library(stringr)
	params <- as.matrix(read.table(filename, sep = "", stringsAsFactors = FALSE, row.names = 1, fill=TRUE))
	param.names <- str_trim(rownames(params))
	
	get.param <- function(param.name){
		param.locale <- which(param.names == param.name)
		if(length(param.locale) == 0){
			return(NA)
			}else{
			param <- params[param.locale,]
			param <- param[which(param != "")]
			final.param <- paste(str_trim(param), collapse = " ")
			return(final.param)
			}
		}
	
	#================================================
	# general parameters
	#================================================
	gen.param <- c("traits", "covariates", "marker.covariates", "traits.scaled", "traits.normalized", "scan.what", "eig.which", "use.kinship", "kinship.type", "locus", "pop", "pval.correction")
	gen.vals <- unlist(lapply(gen.param, get.param))
	
	#================================================
	#single scan parameters
	#================================================
	single.param <- c("ref.allele", "singlescan.perm")
	single.vals <- unlist(lapply(single.param, get.param))

	#================================================
	# marker selection
	#================================================
	marker.param <- c("marker.selection.method", "SNPfile", "peak.density", "tolerance", "window.size", "num.alleles.in.pairscan", "bp.buffer", "organism")
	marker.vals <- unlist(lapply(marker.param, get.param))

	#================================================
	# pair scan
	#================================================
	pair.param <- c("max.pair.cor", "min.per.geno", "pairscan.null.size")
	pair.vals <- unlist(lapply(pair.param, get.param))

	
	param.table <- matrix(c(gen.vals, single.vals, marker.vals, pair.vals), ncol = 1)
	rownames(param.table) <- c(gen.param, single.param, marker.param, pair.param)
	return(param.table)

}
