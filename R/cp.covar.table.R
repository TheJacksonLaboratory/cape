#This script takes a variable from the phenotype matrix
#for example, diet treament or sex and creates a marker
#variable that can be used as a covariate.
#It creates a marker that is numeric and assigns the 
#numeric value to each of the allels at all loci for 
#the given individual.

###############################################################################
#' Method getCovarTable
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod getCovarTable
setGeneric("getCovarTable", function(x) standardGeneric("getCovarTable"))

#' @rdname Cape-class
#' @aliases getCovarTable,Cape-method
setMethod("getCovarTable", "Cape", function(x) x@covar_table)

setGeneric("cpCreateCovarTable<-", function(x, value) standardGeneric("cpCreateCovarTable<-"))

setMethod("cpCreateCovarTable<-", "Cape", function(x, value) {

  # TODO add in this function
	# check.underscore(data.obj)
	
	marker.locale <- get.col.num(x@pheno, value)

	#make a separate covariate table, then modify the dimnames
	#in the genotype object to include the covariates
	#do not modify the genotype object
	
	covar.table <- x@pheno[,marker.locale,drop=FALSE]
	rownames(covar.table) <- rownames(x@pheno)
	x@covar_table <- covar.table
	
		
	#take the phenotypes made into markers out of the phenotype matrix
	x@pheno <- x@pheno[,-marker.locale]
	x@non_allelic_covar <- value
	x@geno_names[[3]] <- c(x@geno_names[[3]], value)
	x@chromosome <- c(x@chromosome, rep(0, length(value)))
	x@marker_location <- c(x@marker_location, 1:length(value))
	
	x
	
})