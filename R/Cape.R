#' The CAPE data object
#' 
#' Class \code{Cape} defines a CAPE analysis object.
#'
#' @name Cape-class
#' @rdname Cape-class
#' @exportClass Cape
#'
#' @slot parameter_file string, full path to YAML file with initialization parameters
#' @slot results_path string, full path to directory for storing results (optional, a directory will be created if one is not specified)
#' @slot pheno A phenotype matrix
#' @slot chromosome A chromosome character list
#' @slot marker_num An integer list of marker numbers along a chromosome
#' @slot marker_location A numeric list of positions in centiMorgans
#' @slot geno_names A list of character names for each genotype, e.g., c("A", "B")
#' @slot geno An array where the dimension names must be "sample", "allele", and "locus"
#' @slot ref_allele A character from the geno_names that represents the wild type
#' @slot covar_table a matrix of 
#' @slot flat_geno a flattened genotype matrix
#' @slot non_allelic_covar covariate
#' @slot p_covar
#' @slot g_covar
#' @slot p_covar_table
#' @slot g_covar_table
#' @slot model_family
#' @slot scan_what string, "eigentraits", "normalized.traits" or "raw.traits"
#' @slot ET eigentraits
#' @slot singular_values
#' @slot right_singular_vectors
#' @slot traits_scaled boolean
#' @slot traits_normalized boolean
#' 
#' @export
Cape <- R6::R6Class(
  "Cape",
  portable = FALSE,
  class = FALSE,
  cloneable = FALSE,
  public = list(
    parameter_file = NULL,
    results_path = NULL,
    pheno = NULL,
    chromosome = NULL,
    marker_num = NULL,
    marker_location = NULL,
    geno_names = NULL,
    geno = NULL,
    ref_allele = NULL,
    covar_table = NULL,
    flat_geno = NULL,
    non_allelic_covar = NULL,
    p_covar = NULL,
    g_covar = NULL,
    p_covar_table = NULL,
    g_covar_table = NULL,
    model_family = NULL,
    scan_what = NULL,
    ET = NULL,
    singular_values = NULL,
    right_singular_vectors = NULL,
    traits_scaled = NULL,
    traits_normalized = NULL,
    
    # this function assigns variables from the parameter file
    # to attributes in the Cape object
    assign_parameters = function() {
      
      parameter.table <- read.parameters(self$parameter_file)
      for(name in names(parameter.table)){
        val <- parameter.table[[name]]
        self[[name]] <- val
      }
    },
    initialize = function(
      parameter_file = NULL,
      results_path = NULL,
      pheno = NULL,
      chromosome = NULL,
      marker_num = NULL,
      marker_location = NULL,
      geno_names = NULL,
      geno = NULL,
      ref_allele = NULL,
      covar_table = NULL,
      flat_geno = NULL,
      non_allelic_covar = NULL,
      p_covar = NULL,
      g_covar = NULL,
      p_covar_table = NULL,
      g_covar_table = NULL,
      model_family = NULL,
      scan_what = NULL,
      ET = NULL,
      singular_values = NULL,
      right_singular_vectors = NULL,
      traits_scaled = NULL,
      traits_normalized = NULL
    ) {
      self$parameter_file <- parameter_file
      if (missing(results_path)) {
        # if the path isn't suplied, take the parameter file's name and append
        # the date and time to create the results directory
        param_name <- tools::file_path_sans_ext(basename(self$parameter_file))
        dt <- format(Sys.time(), "%Y%m%d_%H%M")
        results_path <- paste(param_name, dt, sep = "_")
        self$results_path <- file.path(".", results_path)
        
      } else {
        self$results_path <- results_path
      }
      dir.create(self$results_path, showWarnings = FALSE)
      self$results_path <- results_path
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
      self$covar_table <- covar_table
      self$flat_geno <- flat_geno
      self$non_allelic_covar <- non_allelic_covar
      self$p_covar <- p_covar
      self$g_covar <- g_covar
      self$p_covar_table <- p_covar_table
      self$g_covar_table <- g_covar_table
      self$model_family <- model_family
      self$scan_what <- scan_what
      self$ET <- ET
      self$singular_values <- singular_values
      self$right_singular_vectors <- right_singular_vectors
      self$traits_scaled <- traits_scaled
      self$traits_normalized <- traits_normalized
      # assign parameters from the parameter_file
      self$assign_parameters()
      check.underscore(self)
      check.bad.markers(self)
    },
    plot_svd = function(filename) {
      
      svd.file <- file.path(self$results_path, filename)
      switch(
        tolower(tools::file_ext(filename)),
        "pdf" = pdf(svd.file, width = 7, height = 7),
        "png" = png(svd.file, res = 300, width = 7, height = 7, units = "in"),
        "jpeg" = jpeg(svd.file, res = 300, width = 7, height = 7, units = "in"),
        "jpg" = jpeg(svd.file, res = 300, width = 7, height = 7, units = "in")
      )
      plotSVD(data.obj, orientation = "vertical", show.var.accounted = TRUE)
      dev.off()
      
    },
    set_pheno = function(val) {
      self$pheno <- val
      invisible(self)
    },
    set_geno = function(val) {
      self$geno <- val
      invisible(self)
    },
    create_covar_table = function(value) {
      
      marker.locale <- get.col.num(self$pheno, value)
      
      # make a separate covariate table, then modify the dimnames
      # in the genotype object to include the covariates
      # do not modify the genotype object
      
      self$covar_table <- self$pheno[,marker.locale,drop=FALSE]
      rownames(self$covar_table) <- rownames(self$pheno)
      
      
      # take the phenotypes made into markers out of the phenotype matrix
      self$pheno <- self$pheno[,-marker.locale]
      self$non_allelic_covar <- value
      self$geno_names[[3]] <- c(self$geno_names[[3]], value)
      self$chromosome <- c(self$chromosome, rep(0, length(value)))
      self$marker_location <- c(self$marker_location, 1:length(value))
      invisible(self)
    }
  ),
  lock_objects = FALSE,
  lock_class = TRUE
)
