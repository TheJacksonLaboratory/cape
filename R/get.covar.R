#' Get the phenotype-genotype covariate table
#' 
#' This function returns a covariate table that combines
#' the covariates specified from phenotypes and those
#' specified from genotypes
#'
#' @param data.obj a \code{\link{Cape}} object
#'
#' @return Returns a list with "covar.names", "covar.type", "covar.loc", "covar.table"
#' 
#' @export
get.covar <- function(data.obj){
		
	covar.table <- cbind(data.obj$p_covar_table, data.obj$g_covar_table)
	
	if(!is.null(data.obj$p_covar_table)){
		num.p.covar <- dim(data.obj$p_covar_table)
		if (length(num.p.covar) == 1) {
		  # there's only one covar column
		  num.p.covar <- 1
		} else {
		  num.p.covar <- num.p.covar[2]
		}
		}else{
		num.p.covar <- 0	
		}
		
	if(!is.null(data.obj$g_covar_table)){
		num.g.covar <- dim(data.obj$g_covar_table)[2]
		}else{
		num.g.covar <- 0	
		}
	
	covar.type <- c(rep("p", num.p.covar), rep("g", num.g.covar))

	covar <- as.vector(c(data.obj$p_covar, data.obj$g_covar[1,]))
	
		covar.names <- c(data.obj$p_covar, data.obj$g_covar[1,])
		p.covar.loc <- NULL
		g.covar.loc <- NULL
				
		p.covar.locale <- which(covar.type == "p")
		g.covar.locale <- which(covar.type == "g")
		
		if(length(p.covar.locale) > 0){
			p.covar.loc <- rep(1:length(which(covar.type == "p")))
			}
		if(length(g.covar.locale) > 0){
			g.covar.loc <- data.obj$g_covar["position",]
			}
		covar.loc <- c(p.covar.loc, g.covar.loc)
		covar.table <- cbind(data.obj$p_covar_table, data.obj$g_covar_table)
		
		if(length(covar.loc) < length(covar)){
			not.found <- setdiff(covar, covar.names)
			if(length(not.found) > 0){
				cat("I could not find the following covariates:\n")
				cat(not.found, sep = "\n")
				stop()
				}
			}
	
	result <- list("covar.names" = covar.names, "covar.type" = covar.type, "covar.loc" = covar.loc, "covar.table" = covar.table)
	return(result)	
	
}