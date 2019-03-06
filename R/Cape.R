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
#' @slot geno An array where the dimension names must be "sample", "allele", and "locus"
#' @slot ref_allele A character from the geno_names that represents the wild type
#' @slot parameters TODO mayhap we should change this?
#' @slot covar_table a matrix of 
#' @slot flat_geno a flattened genotype matrix
#' @slot non_allelic_covar covariate
#' 
#' @export
Cape <- R6::R6Class(
  "Cape",
  public = list(
    pheno = NULL,
    chromosome = NULL,
    marker_num = NULL,
    marker_location = NULL,
    geno_names = NULL,
    geno = NULL,
    ref_allele = NULL,
    parameters = NULL,
    covar_table = NULL,
    flat_geno = NULL,
    non_allelic_covar = NULL,
    
    initialize = function(pheno = NA, chromosome = NA, marker_num = NA,
                          marker_location = NA, geno = NA, geno_names = NA,
                          parameters = NA, ref_allele = NA, covar_table = NA,
                          flat_geno = NA, non_allelic_covar = NA) {
      self$pheno <- pheno
      self$chromosome <- chromosome
      self$marker_num <- marker_num
      self$marker_location <- marker_location
      self$geno <- geno
      self$geno_names <- geno_names
      if (!missing(ref_allele)) {
        stopifnot(is.character(ref_allele))
        self$ref_allele <- ref_allele
      }
      self$parameters <- parameters
      self$covar_table <- covar_table
      self$flat_geno <- flat_geno
      self$non_allelic_covar <- non_allelic_covar
    },
    set_pheno = function(val) {
      self$pheno <- val
    },
    set_geno = function(val) {
      self$geno <- val
    },
    create_covar_table = function(value) {
      marker.locale <- get.col.num(self$pheno, value)
      
      #make a separate covariate table, then modify the dimnames
      #in the genotype object to include the covariates
      #do not modify the genotype object
      
      covar.table <- self$pheno[,marker.locale,drop=FALSE]
      rownames(covar.table) <- rownames(self$pheno)
      self$covar_table <- covar.table
      
      
      #take the phenotypes made into markers out of the phenotype matrix
      self$pheno <- self$pheno[,-marker.locale]
      self$non_allelic_covar <- value
      self$geno_names[[3]] <- c(self$geno_names[[3]], value)
      self$chromosome <- c(self$chromosome, rep(0, length(value)))
      self$marker_location <- c(self$marker_location, 1:length(value))
      invisible(self)
    }
  )
)
