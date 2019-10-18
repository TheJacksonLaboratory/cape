#' The CAPE data object
#'
#' Class \code{Cape} defines a CAPE analysis object.
#'
#' @name Cape-class
#' @rdname Cape-class
#' @exportClass Cape
#'
#' @slot parameter_file string, full path to YAML file with initialization
#'   parameters
#' @slot results_path string, full path to directory for storing results
#'   (optional, a directory will be created if one is not specified)
#' @slot save_results boolean, default: TRUE
#' @slot use_saved_results boolean, default: TRUE
#' @slot orginism options are "mouse" or "human"
#' @slot pheno A phenotype matrix
#' @slot chromosome A chromosome character list
#' @slot marker_num An integer list of marker numbers along a chromosome
#' @slot marker_location A numeric list of positions in centiMorgans
#' @slot marker_selection_method options are "top.effects", "from.list",
#'   "uniform", "by.gene", "effects.dist"
#' @slot bp_buffer when the marker selection method is "by.gene" this finds genes 
#'   from Ensemble that are within this number of base pairs of the input marker's 
#'   position
#' @slot geno_names A list of character names for each genotype, e.g., c("A","B")
#' @slot geno An array where the dimension names must be "sample", "allele", and "locus"
#' @slot geno_for_pairscan
#' @slot effect_size_cutoff from select.markers.for.pairscan.R
#' @slot peak_density from select.markers.for.pariscan.R
#' @slot window_size from select.markers.for.pariscan.R
#' @slot tolerance from select.markers.for.pairscan.R
#' @slot ref_allele A character from the geno_names that represents the wildtype
#' @slot covar_table a matrix of
#' @slot flat_geno a flattened genotype matrix
#' @slot non_allelic_covar covariate
#' @slot num_alleles_in_pairscan
#' @slot max_pair_cor the maximum Pearson correlation that two markers are allowed
#' @slot min_per_genotype minimum number of individuals allowable per genotype
#' @slot pairscan_null_size
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
#' @slot var_to_var_influences_perm added in \code{\link{error.prop}} 
#' if error propogation is performed on permuted test statistics
#' @slot var_to_var_influences added in \code{\link{error.prop}} 
#' if error propogation is performed on un-permuted test statistics
#' @slot pval_correction options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
#' @slot linkage_blocks_collapsed see \code{\link{get.network}}
#' @slot linkage_blocks_full see \code{\link{get.network}}
#' @slot var_to_var_p_val see \code{\link{get.network}}
#' @slot max_var_to_pheno_influence see \code{\link{get.network}}
#' @slot collapsed_net see \code{\link{get.network}}
#' @slot full_net see \code{\link{get.network}}
#' @slot use_kinship boolean
#'
#' @export
Cape <- R6::R6Class(
  "Cape",
  portable = FALSE,
  class = FALSE,
  cloneable = FALSE,
  private = list(
    .geno_for_pairscan = NULL,
    .marker_selection_method = NULL,
    .linkage_blocks_collapsed = NULL,
    .collapsed_net = NULL
  ),
  active = list(
    geno_for_pairscan = function(value) {
      if (missing(value)) {
        private$.geno_for_pairscan
      } else {
        private$.geno_for_pairscan <- value
        self
      }
    },
    marker_selection_method = function(value) {
      if (missing(value)) {
        private$.marker_selection_method
      } else {
        private$.marker_selection_method <- value
        self
      }
    },
    linkage_blocks_collapsed = function(value) {
      if (missing(value)) {
        private$.linkage_blocks_collapsed
      } else {
        private$.linkage_blocks_collapsed <- value
        self
      }
    },
    linkage_blocks_full = function(value) {
      if (missing(value)) {
        private$.linkage_blocks_full
      } else {
        private$.linkage_blocks_full <- value
        self
      }
    },
    collapsed_net = function(value) {
      if (missing(value)) {
        private$.collapsed_net
      } else {
        private$.collapsed_net <- value
        self
      }
    }
  ),
  public = list(
    parameter_file = NULL,
    results_path = NULL,
    save_results = NULL,
    use_saved_results = NULL,
    organism = NULL,
    pheno = NULL,
    chromosome = NULL,
    marker_num = NULL,
    marker_location = NULL,
    bp_buffer = NULL,
    geno_names = NULL,
    geno = NULL,
    effect_size_cutoff = NULL,
    peak_density = NULL,
    window_size = NULL,
    tolerance = NULL,
    ref_allele = NULL,
    covar_table = NULL,
    flat_geno = NULL,
    non_allelic_covar = NULL,
    num_alleles_in_pairscan = NULL,
    max_pair_cor = NULL,
    min_per_genotype = NULL,
    pairscan_null_size = NULL,
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
    var_to_var_influences_perm = NULL,
    var_to_var_influences = NULL,
    pval_correction = NULL,
    var_to_var_p_val = NULL,
    max_var_to_pheno_influence = NULL,
    full_net = NULL,
    use_kinship = NULL,
    
    # this function assigns variables from the parameter file
    # to attributes in the Cape object
    assign_parameters = function() {
      
      parameter.table <- read.parameters(self$parameter_file)
      for(name in names(parameter.table)){
        val <- parameter.table[[name]]
        self[[name]] <- val
      }
    },
    check_inputs = function() {
      stopifnot(length(self$chromosome) == length(self$marker_location))
      stopifnot(length(self$chromosome) == length(self$marker_num))
      stopifnot(length(self$chromosome) == length(self$geno_names$locus))
    },
    check_geno_names = function() {
      # TODO make sure that individual names match between the pheno object, geno object, and geno names
      stopifnot(TRUE)
      
    },
    initialize = function(
      parameter_file = NULL,
      results_path = NULL,
      save_results = TRUE,
      use_saved_results = TRUE,
      organism = NULL,
      pheno = NULL,
      chromosome = NULL,
      marker_num = NULL,
      marker_location = NULL,
      bp_buffer = NULL,
      geno_names = NULL,
      geno = NULL,
      .geno_for_pairscan = NULL,
      effect_size_cutoff = NULL,
      peak_density = NULL,
      window_size = NULL,
      tolerance = NULL, 
      ref_allele = NULL,
      covar_table = NULL,
      flat_geno = NULL,
      non_allelic_covar = NULL,
      num_alleles_in_pairscan = NULL,
      max_pair_cor = NULL,
      min_per_genotype = NULL,
      pairscan_null_size = NULL,
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
      var_to_var_influences_perm = NULL,
      var_to_var_influences = NULL,
      pval_correction = NULL,
      var_to_var_p_val = NULL,
      max_var_to_pheno_influence = NULL,
      full_net = NULL,
      use_kinship = NULL
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
      self$save_results <- save_results
      self$use_saved_results <- use_saved_results
      self$organism <- organism
      self$pheno <- pheno
      self$chromosome <- chromosome
      self$marker_num <- marker_num
      self$marker_location <- marker_location
      if (missing(bp_buffer)) {
        self$bp_buffer <- 1000
      }
      self$geno <- geno
      self$effect_size_cutoff <- effect_size_cutoff
      self$peak_density <- peak_density
      self$window_size <- window_size
      self$tolerance <- tolerance
      self$geno_names <- geno_names
      if (!missing(ref_allele)) {
        stopifnot(is.character(ref_allele))
        self$ref_allele <- ref_allele
      }
      self$covar_table <- covar_table
      self$flat_geno <- flat_geno
      self$non_allelic_covar <- non_allelic_covar
      self$num_alleles_in_pairscan <- num_alleles_in_pairscan
      self$max_pair_cor <- max_pair_cor
      self$min_per_genotype <- min_per_genotype
      self$pairscan_null_size <- pairscan_null_size
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
      self$var_to_var_influences_perm <- var_to_var_influences_perm
      self$var_to_var_influences <- var_to_var_influences
      self$pval_correction <- pval_correction
      self$var_to_var_p_val <- var_to_var_p_val
      self$max_var_to_pheno_influence <- max_var_to_pheno_influence
      self$full_net <- full_net
      self$use_kinship <- use_kinship
      # assign parameters from the parameter_file
      self$assign_parameters()
      self$check_inputs()
      self$check_geno_names()
      delete.underscore(self)
      check.underscore(self)
      check.bad.markers(self)
      # TODO make sure that individual names match between the pheno object, geno object, and geno names
    },
    plot_svd = function(filename) {
      
      full.path <- file.path(self$results_path, filename)
      
      switch(
        tolower(tools::file_ext(filename)),
        "pdf" = pdf(full.path, width = 7, height = 7),
        "png" = png(full.path, res = 300, width = 7, height = 7, units = "in"),
        "jpeg" = jpeg(full.path, res = 300, width = 7, height = 7, units = "in"),
        "jpg" = jpeg(full.path, res = 300, width = 7, height = 7, units = "in")
      )
      plotSVD(self, orientation = "vertical", show.var.accounted = TRUE)
      dev.off()
      
    },
    plot_singlescan = function(filename, singlescan.obj, width = 20, height = 6, units = "in", res = 300, 
                               standardized = TRUE, allele.labels = NULL, alpha = c(0.05, 0.01), include.covars = TRUE, 
                               line.type = "l", pch = 16, cex = 0.5, lwd = 3, traits = NULL) {
      
      full.path <- file.path(self$results_path, filename)
      
      jpeg(full.path, width = width, height = height, units = units, res = res)
      plotSinglescan(self, singlescan.obj = singlescan.obj, standardized = standardized, allele.labels = allele.labels, 
                     alpha = alpha, include.covars = include.covars, line.type = line.type, pch = pch, cex = cex, 
                     lwd = lwd, traits = traits)
      dev.off()
      
    },
    plot_pairscan = function(filename, pairscan.obj, phenotype = NULL, 
                             show.marker.labels = TRUE, show.alleles = FALSE) {
      
      # filename is usually "Pairscan.Regression.pdf"
      
      full.path <- file.path(self$results_path, filename)
      
      plotPairscan(self, pairscan.obj, phenotype = phenotype, pdf.label = full.path, 
                   show.marker.labels = show.marker.labels, show.alleles = show.alleles)
      
    },
    plot_variant_influences = function(filename, width = 10, height = 7,
                                       p.or.q = p.or.q, standardize = FALSE, 
                                       not.tested.col = "lightgray", 
                                       covar.width = 30, pheno.width = 30) {
      
      full.path <- file.path(self$results_path, filename)
      
      pdf(full.path, width = width, height = height)
      plotVariantInfluences(self, p.or.q = p.or.q, standardize = FALSE, 
                            not.tested.col = "lightgray", 
                            covar.width = 30, pheno.width = 30)
      dev.off()
      
    },
    plot_network_do = function(filename, label.gap = 10, label.cex = 1.5, show.alleles = FALSE) {
      
      full.path <- file.path(self$results_path, filename)
      
      pdf(full.path)
      plotNetworkDO(self, label.gap = label.gap, label.cex = label.cex, show.alleles = show.alleles)
      dev.off()
    },
    plot_full_network = function(filename, zoom = 1.2, node.radius = 0.3, label.nodes = TRUE, label.offset = 0.4, label.cex = 0.5, 
                                 bg.col = "lightgray", arrow.length = 0.1, layout.matrix = "layout_with_kk", legend.position = "topright", 
                                 edge.lwd = 1, legend.radius = 2, legend.cex = 0.7, xshift = -1) {
      
      full.path <- file.path(self$results_path, filename)
      
      pdf(full.path)
      plotFullNetwork(self, zoom = zoom, node.radius = node.radius, label.nodes = label.nodes, label.offset = label.offset, label.cex = label.cex, 
                      bg.col = bg.col, arrow.length = arrow.length, layout.matrix = layout.matrix, legend.position = legend.position, 
                      edge.lwd = edge.lwd, legend.radius = legend.radius, legend.cex = legend.cex, xshift = xshift)
      dev.off()
    },
    write_variant_influences = function(filename, p.or.q = 0.05) {
      
      full.path <- file.path(self$results_path, filename)
      
      writeVariantInfluences(self, p.or.q = max(c(p.or.q, 0.2)), filename = full.path)
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
    },
    save_rds = function(object, filename) {
      # only save the results RData file if save_results is TRUE
      if (self$save_results) {
        full.path <- file.path(self$results_path, filename)
        saveRDS(object, full.path)
      }
    },
    read_rds = function(filename) {
      full.path <- file.path(self$results_path, filename)
      # only return the results RData file if use_saved_results is TRUE
      if ((self$use_saved_results) & (file.exists(full.path))) {
        return(readRDS(full.path))
      } else {
        return(FALSE)
      }
    }
  ),
  lock_objects = FALSE,
  lock_class = TRUE
)
