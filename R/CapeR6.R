CapeR6 <- R6Class(
  "CapeR6",
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
    
    initialize = function(pheno = NA, geno = NA, covar_table = NA) {
      self$pheno <- pheno
      self$geno <- geno
    },
    set_pheno = function(val) {
      self$pheno <- val
    },
    set_geno = function(val) {
      self$geno <- val
    },
    cteate_covar_table = function(value) {
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
