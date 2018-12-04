#' The CAPE data object
#' 
#' Class \code{Cape} defines a CAPE analysis object.
#'
#' @name Cape-class
#' @rdname Cape-class
#' @exportClass Cape
#'
#' @slot pheno A phenotype matrix
#' @slot chromosome A chromosome character list
#' @slot marker_num An integer list of marker numbers along a chromosome
#' @slot marker_location A numeric list of positions in centiMorgans
#' @slot geno_names A list of character names for each genotype, e.g., c("A", "B")
#' @slot ref_allele A character from the geno_names that represents the wild type
#' @slot geno An array where the dimension names must be "sample", "allele", and "locus"
#' @slot covar.table a matrix of 
#' @slot parameters TODO mayhap we should change this?
setClass("Cape", 
  slots = c(
    pheno = "matrix", 
    chromosome = "character",
    marker_num = "integer",
    marker_location = "numeric",
    geno_names = "list",
    ref_allele = "character",
    geno = "array",
    parameters = "array"
  )
)

setValidity("Cape",
    function(object)
    {
      cl <- class(object)
      for(i in slotNames(cl)) {
        if ( length(slot(object,i)) < 1 ) {
          return(paste("The ", i, " slot cannot be null."))
        }
      }
      
      # TODO add a validation check on the geno object's dimension names
      # TODO they should be "sample", "allele", and "locus"
      
      # TODO add a validation check on the genotype names
      # gn <- object@geno_names$sample
      #   gs <- dimnames(object@geno)$sample
      #   if ( length(intersect(gs, gn)) <= 1 )
      #       return("The sample names in @geno and @geno_names don't match.")
        TRUE
    }
)

###############################################################################
#' Method getPheno
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod getPheno
setGeneric("getPheno", function(x) standardGeneric("getPheno"))

#' @rdname Cape-class
#' @aliases getPheno,Cape-method
setMethod("getPheno", "Cape", function(x) x@pheno)

###############################################################################
#' Method setPheno
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod setPheno<-
setGeneric("setPheno<-", function(x, value) standardGeneric("setPheno<-"))

#' @rdname Cape-class
#' @aliases setPheno,Cape-method
#' @param value a phenotype matrix
setMethod("setPheno<-", "Cape", function(x, value) {
    x@pheno <- value
    x
})

###############################################################################
#' Method getChromosome
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod getChromosome
setGeneric("getChromosome", function(x) standardGeneric("getChromosome"))

#' @rdname Cape-class
#' @aliases getChromosome,Cape-method
setMethod("getChromosome", "Cape", function(x) x@chromosome)

###############################################################################
#' Method setChromosome
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod setChromosome<-
setGeneric("setChromosome<-", function(x, value) standardGeneric("setChromosome<-"))

#' @rdname Cape-class
#' @aliases setChromosome,Cape-method
#' @param value a chromosome character list
setMethod("setChromosome<-", "Cape", function(x, value) {
    x@chromosome <- value
    x
})

###############################################################################
#' Method get
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod getMarkerNum
setGeneric("getMarkerNum", function(x) standardGeneric("getMarkerNum"))

#' @rdname Cape-class
#' @aliases getMarkerNum,Cape-method
setMethod("getMarkerNum", "Cape", function(x) x@marker_num)

###############################################################################
#' Method setMarkerNum
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod setMarkerNum<-
setGeneric("setMarkerNum<-", function(x, value) standardGeneric("setMarkerNum<-"))

#' @rdname Cape-class
#' @aliases setMarkerNum,Cape-method
#' @param value An integer list of marker numbers along a chromosome
setMethod("setMarkerNum<-", "Cape", function(x, value) {
    x@marker_num <- value
    x
})

###############################################################################
#' Method getMarkerLocation
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod getMarkerLocation
setGeneric("getMarkerLocation", function(x) standardGeneric("getMarkerLocation"))

#' @rdname Cape-class
#' @aliases getMarkerLocation,Cape-method
setMethod("getMarkerLocation", "Cape", function(x) x@marker_location)

###############################################################################
#' Method setCMarkerLocation
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod setMarkerLocation<-
setGeneric("setMarkerLocation<-", function(x, value) standardGeneric("setMarkerLocation<-"))

#' @rdname Cape-class
#' @aliases setMarkerLocation,Cape-method
#' @param value A numeric list of positions in centiMorgans
setMethod("setMarkerLocation<-", "Cape", function(x, value) {
    x@marker_location <- value
    x
})

###############################################################################
#' Method getGenoNames
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod getGenoNames
setGeneric("getGenoNames", function(x) standardGeneric("getGenoNames"))

#' @rdname Cape-class
#' @aliases getGenoNames,Cape-method
setMethod("getGenoNames", "Cape", function(x) x@geno_names)

###############################################################################
#' Method setGenoNames
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod setGenoNames<-
setGeneric("setGenoNames<-", function(x, value) standardGeneric("setGenoNames<-"))

#' @rdname Cape-class
#' @aliases setGenoNames,Cape-method
#' @param value A list of character names for each genotype, e.g., c("A", "B")
setMethod("setGenoNames<-", "Cape", function(x, value) {
    x@geno_names <- value
    x
})

###############################################################################
#' Method getRefAllele
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod getRefAllele
setGeneric("getRefAllele", function(x) standardGeneric("getRefAllele"))

#' @rdname Cape-class
#' @aliases getRefAllele,Cape-method
setMethod("getRefAllele", "Cape", function(x) x@ref_allele)

###############################################################################
#' Method setRefAllele
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod setRefAllele<-
setGeneric("setRefAllele<-", function(x, value) standardGeneric("setRefAllele<-"))

#' @rdname Cape-class
#' @aliases set,Cape-method
#' @param value A character from the geno_names that represents the wild type
setMethod("setRefAllele<-", "Cape", function(x, value) {
    x@ref_allele <- value
    x
})

###############################################################################
#' Method getGeno
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod getGeno
setGeneric("getGeno", function(x) standardGeneric("getGeno"))

#' @rdname Cape-class
#' @aliases getGeno,Cape-method
setMethod("getGeno", "Cape", function(x) x@geno)

###############################################################################
#' Method setGeno
#' @name Cape-class
#' @rdname Cape-class
#' @exportMethod setGeno<-
setGeneric("setGeno<-", function(x, value) standardGeneric("setGeno<-"))

#' @rdname Cape-class
#' @aliases setGeno,Cape-method
#' @param value An array where the dimension names must be "sample", "allele", and "locus"
setMethod("setGeno<-", "Cape", function(x, value) {
    x@geno <- value
    x
})
