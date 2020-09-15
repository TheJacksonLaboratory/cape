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
#' @slot yaml_parameters string representing YAML CAPE parameters. See the
#'   vignette for more descriptions of individual parameter settings.
#' @slot results_path string, full path to directory for storing results
#'   (optional, a directory will be created if one is not specified)
#' @slot save_results Whether to save cape results. Defaults to TRUE.
#' @slot use_saved_results Whether to use existing results from a 
#'   previous run. This can save time if re-running an analysis, but
#'   can lead to problems if the old run and new run have competing settings.
#'   If errors arise, and use_saved_results is set to TRUE, try setting it 
#'   to FALSE, or deleting previous results.
#' @slot pheno A matrix containing the traits to be analyzed. Traits are in
#'    columns and individuals are in rows.
#' @slot chromosome A vector the same length as the number of markers indicating
#'    which chromosome each marker lives on.
#' @slot marker_num A vector the same length as the number of markers indicating
#'    the index of each marker
#' @slot marker_location A vector the same length as the number of markers indicating
#'    the genomic position of each marker. The positions are primarily used for plotting
#'    and can be in base pairs, centiMorgans, or dummy variables.
#' @slot marker_selection_method A string indicating how markers should be 
#'   selected for the pairscan. Options are "top.effects" or "from.list."
#'   If "top.effects," markers are selected using main effect sizes. 
#'   If "from.list," markers are specified using a vector of marker names. 
#'   See \code{\link{select.markers.for.pairscan}}.
#' @slot geno_names The dimnames of the genotype array. The genotype array is a three-dimensional
#'   array in which rows are individuals, columns are alleles, and the third dimension houses
#'   the markers. Genotypes are pulled for analysis using \code{\link{get.geno}} based on
#'   geno_names. Only the individuals, alleles, and markers listed in geno_names are
#'   taken from the genotype matrix. Functions that remove markers and individuals from
#'   analysis always operate on geno_names in addition to other relevant slots.
#'   The names of geno_names must be "mouse", "allele", "locus."
#' @slot geno A three dimentional array holding genotypes for each animal for each allele
#'   at each marker. The genotypes are continuously valued probabilities ranging from 0 to 1. 
#'   The dimnames of geno must be "mouse", "allele", and "locus," even if the individuals are
#'   not mice.
#' @slot geno_for_pairscan A two-dimensional matrix holding the genotypes that will be analyzed
#'   in the pairscan. Alleles are in columns and individuals are in rows. As in the geno array, 
#'   values are continuous probabilities ranging from 0 to 1.
#' @slot peak_density The density parameter for \code{\link{select.markers.for.pairscan}}.
#'   Determines how densely markers under an individual effect size peak are selected 
#'   for the pairscan if marker_selection_method is TRUE. Defaults to 0.5.
#' @slot window_size The window size used by \code{\link{select.markers.for.pairscan}}.
#'   It specifies how many markers are used to smooth effect size curves for automatic peak
#'   identification. If set to NULL, window_size is determined automatically. Used when 
#'   marker_selection_method is TRUE.
#' @slot tolerance The wiggle room afforded to \code{\link{select.markers.for.pairscan}} in 
#'   finding a target number of markers. If num_alleles_in_pairscan is 100 and the tolerance 
#'   is 5, the algorithm will stop when it identifies anywhere between 95 and 105 markers 
#'   for the pairscan.
#' @slot ref_allele A string of length 1 indicating which allele to use as the reference allele.
#'   In two-parent crosses, this is usually allele A. In DO/CC populations, we recommend using
#'   B as the reference allele. B is the allele from the C57Bl6/J mouse, which is often used as
#'   a reference strain.
#' @slot alpha The significance level for calculating effect size thresholds in the 
#'   \code{\link{singlescan}}. If singlescan_perm is 0, this parameter is ignored.
#' @slot covar_table A matrix of covariates with covariates in columns and individuals
#'   in rows. Must be numeric.
#' @slot num_alleles_in_pairscan The number of alleles to test in the pairwise scan. 
#'   Because Cape is computationally intensive, we usually need to test only a subset
#'   of available markers in the pairscan, particularly if the kinship correction is
#'   being used.
#' @slot max_pair_cor the maximum Pearson correlation between two markers. If their
#'   correlation exceeds this value, they will not be tested against each other in the
#'   pairscan. This threshold is set to prevent false positive arising from testing
#'   highly correlated markers. If this value is set to NULL, min_per_genotype must
#'   be specified.
#' @slot min_per_genotype minimum The minimum number of individuals allowable per
#'   genotype combination in the pair scan. If for a given marker pair, one of the 
#'   genotype combinations is underrepresented, the marker pair is not tested. If 
#'   this value is NULL, max_pair_cor must be specified.
#' @slot pairscan_null_size The total size of the null distribution.
#'   This is DIFFERENT than the number of permutations to run. Each permutation
#'   generates n choose 2 elements for the pairscan. So for example, a permutation
#'   that tests 100 pairs of markers will generate a null distribution of size 4950.
#'   This process is repeated until the total null size is reached. If the null size
#'   is set to 5000, two permutations of 100 markers would be done to get to a null
#'   distribution size of 5000.
#' @slot p_covar A vector of strings specifying the names of covariates derived
#'   from traits. See \code{\link{pheno2covar}}.
#' @slot g_covar A vector of strings specifying the names of covariates derived 
#'   from genetic markers. See \code{\link{marker2covar}}.
#' @slot p_covar_table A matrix holding the individual values for each
#'   trait-derived covariate. See \code{\link{pheno2covar}}.
#' @slot g_covar_table A matrix holding the individual values for each 
#'   marker-derived covariate. See \code{\link{marker2covar}}.
#' @slot model_family Indicates the model family of the phenotypes
#'   This can be either "gaussian" or "binomial". If this argument
#'   is length 1, all phenotypes will be assigned to the same
#'   family. Phenotypes can be assigned different model families by
#'   providing a vector of the same length as the number of phenotypes,
#'   indicating how each phenotype should be modeled. See \code{\link{singlescan}}.
#' @slot scan_what A string indicating whether "eigentraits", "normalized.traits", or 
#'   "raw.traits" should be analyzed. See \code{\link{get.pheno}}.
#' @slot ET A matrix holding the eigentraits to be analyzed.
#' @slot singular_values Added by \code{\link{get.eigentraits}}. A vector holding 
#'   the singular values from the singular
#'   value decomposition of the trait matrix. They are used in rotating the 
#'   final direct influences back to trait space from eigentrait space. See
#'   \code{\link{get.eigentraits}} and \code{\link{direct.influence}}.
#' @slot right_singular_vectors Added by \code{\link{get.eigentraits}}. A matrix 
#'   containing the right singular vectors from the singular
#'   value decomposition of the trait matrix. They are used in rotating the 
#'   final direct influences back to trait space from eigentrait space. See
#'   \code{\link{get.eigentraits}} and \code{\link{direct.influence}}.
#' @slot traits_scaled Whether the traits should be mean-centered and standardized
#'   before analyzing.
#' @slot traits_normalized Whether the traits should be rank Z normalized before
#'   analyzing.
#' @slot var_to_var_influences_perm added in \code{\link{error.prop}} 
#'  The list of results from the error propagation of permuted coefficients.
#' @slot var_to_var_influences added in \code{\link{error.prop}} 
#'  The list of results from the error propagation of coefficients.
#' @slot pval_correction Options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
#' @slot linkage_blocks_collapsed A list containing assignments of markers to linkage blocks 
#'  calculated by \code{\link{linkage.blocks.network}} and \code{\link{plotNetwork}}.
#'  In this list there can be multiple markers assigned to a single linkage block.
#' @slot linkage_blocks_full A list containing assignments of markers to linkage blocks 
#'  when no linkage blocks are calculated. In this list there can only be one marker
#'  per "linkage block". See \code{\link{linkage.blocks.network}} and \code{\link{plotNetwork}}.
#' @slot var_to_var_p_val The final table of cape interaction results calculated by \code{\link{error.prop}}.
#' @slot max_var_to_pheno_influence The final table of cape direct influences of markers to traits
#'  calculated by \code{\link{direct.influence}}.
#' @slot collapsed_net An adjacency matrix holding significant cape interactions between
#'  linkage blocks. See \code{\link{plotNetwork}} and \code{\link{get.network}}.
#' @slot full_net An adjacency matrix holding significant cape interactions between
#'  individual markers. See \code{\link{plotNetwork}} and \code{\link{get.network}}.
#' @slot use_kinship Whether to use a kinship correction in the analysis.
#' @slot kinship_type Which type of kinship matrix to use. Either "overall" 
#' for the overall kinship matrix or "ltco" for leave-two-chromosomes-out.
#' @slot transform_to_phenospace whether to transform to phenospace or not.
#'
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
    #' @field geno_for_pairscan Set \code{\link{.geno_for_pairscan}}.
    geno_for_pairscan = function(value) {
      if (missing(value)) {
        private$.geno_for_pairscan
      } else {
        private$.geno_for_pairscan <- value
        self
      }
    },
    #' @field marker_selection_method Set \code{\link{.marker_selection_method}}.
    marker_selection_method = function(value) {
      if (missing(value)) {
        private$.marker_selection_method
      } else {
        private$.marker_selection_method <- value
        self
      }
    },
    #' @field linkage_blocks_collapsed Set \code{\link{.linkage_blocks_collapsed}}.
    linkage_blocks_collapsed = function(value) {
      if (missing(value)) {
        private$.linkage_blocks_collapsed
      } else {
        private$.linkage_blocks_collapsed <- value
        self
      }
    },
    #' @field linkage_blocks_full Set \code{\link{.linkage_blocks_full}}.
    linkage_blocks_full = function(value) {
      if (missing(value)) {
        private$.linkage_blocks_full
      } else {
        private$.linkage_blocks_full <- value
        self
      }
    },
    #' @field collapsed_net Set \code{\link{.collapsed_net}}.
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
    #' @field parameter_file full path to YAML file with initialization parameters.
    parameter_file = NULL,
    #' @field yaml_parameters string representing YAML CAPE parameters. See the vignette for more descriptions 
    #' of individual parameter settings.
    yaml_parameters = NULL,
    #' @field results_path string, full path to directory for storing results (optional, a directory will be 
    #' created if one is not specified).
    results_path = NULL,
    #' @field save_results Whether to save cape results. Defaults to TRUE.
    save_results = NULL,
    #' @field use_saved_results Whether to use existing results from a previous run. This can save time if 
    #' re-running an analysis, but can lead to problems if the old run and new run have competing settings.
    #' If errors arise, and use_saved_results is set to TRUE, try setting it to FALSE, or deleting previous results.
    use_saved_results = NULL,
    #' @field pheno A matrix containing the traits to be analyzed. Traits are in columns and individuals are in rows.
    pheno = NULL,
    #' @field chromosome A vector the same length as the number of markers indicating which chromosome each marker lives on.
    chromosome = NULL,
    #' @field marker_num A vector the same length as the number of markers indicating the index of each marker.
    marker_num = NULL,
    #' @field marker_location A vector the same length as the number of markers indicating the genomic position of each marker.
    #' The positions are primarily used for plotting and can be in base pairs, centiMorgans, or dummy variables.
    marker_location = NULL,
    #' @field geno_names The dimnames of the genotype array. The genotype array is a three-dimensional array in which rows 
    #' are individuals, columns are alleles, and the third dimension houses the markers. Genotypes are pulled for analysis 
    #' using \code{\link{get.geno}} based on geno_names. Only the individuals, alleles, and markers listed in geno_names are
    #' taken from the genotype matrix. Functions that remove markers and individuals from analysis always operate on geno_names 
    #' in addition to other relevant slots. The names of geno_names must be "mouse", "allele", "locus."
    geno_names = NULL,
    #' @field geno A three dimentional array holding genotypes for each animal for each allele at each marker. The genotypes
    #' are continuously valued probabilities ranging from 0 to 1. The dimnames of geno must be "mouse", "allele", and "locus,"
    #' even if the individuals are not mice.
    geno = NULL,
    #' @field peak_density The density parameter for \code{\link{select.markers.for.pairscan}}.  Determines how densely
    #' markers under an individual effect size peak are selected for the pairscan if marker_selection_method is TRUE.
    #' Defaults to 0.5.
    peak_density = NULL,
    #' @field window_size The window size used by \code{\link{select.markers.for.pairscan}}. It specifies how many markers
    #' are used to smooth effect size curves for automatic peak identification. If set to NULL, window_size is determined
    #' automatically. Used when marker_selection_method is TRUE.
    window_size = NULL,
    #' @field tolerance The wiggle room afforded to \code{\link{select.markers.for.pairscan}} in finding a target number
    #' of markers. If num_alleles_in_pairscan is 100 and the tolerance is 5, the algorithm will stop when it identifies
    #' anywhere between 95 and 105 markers for the pairscan.
    tolerance = NULL,
    #' @field ref_allele A string of length 1 indicating which allele to use as the reference allele. In two-parent crosses,
    #' this is usually allele A. In DO/CC populations, we recommend using B as the reference allele. B is the allele from
    #' the C57Bl6/J mouse, which is often used as a reference strain.
    ref_allele = NULL,
    #' @field alpha The significance level for calculating effect size thresholds in the \code{\link{singlescan}}.
    #' If singlescan_perm is 0, this parameter is ignored.
    alpha = NULL,
    #' @field covar_table A matrix of covariates with covariates in columns and individuals in rows. Must be numeric.
    covar_table = NULL,
    #' @field num_alleles_in_pairscan The number of alleles to test in the pairwise scan. Because Cape is computationally
    #' intensive, we usually need to test only a subset of available markers in the pairscan, particularly if the kinship
    #' correction is being used.
    num_alleles_in_pairscan = NULL,
    #' @field max_pair_cor The maximum Pearson correlation between two markers. If their correlation exceeds this value,
    #' they will not be tested against each other in the pairscan. This threshold is set to prevent false positive arising
    #' from testing highly correlated markers. If this value is set to NULL, min_per_genotype must be specified.
    max_pair_cor = NULL,
    #' @field min_per_genotype minimum The minimum number of individuals allowable per genotype combination in the pair scan.
    #' If for a given marker pair, one of the genotype combinations is underrepresented, the marker pair is not tested.
    #' If this value is NULL, max_pair_cor must be specified.
    min_per_genotype = NULL,
    #' @field pairscan_null_size The total size of the null distribution. This is DIFFERENT than the number of permutations
    #' to run. Each permutation generates n choose 2 elements for the pairscan. So for example, a permutation that tests 100
    #' pairs of markers will generate a null distribution of size 4950. This process is repeated until the total null size
    #' is reached. If the null size is set to 5000, two permutations of 100 markers would be done to get to a null
    #' distribution size of 5000.
    pairscan_null_size = NULL,
    #' @field p_covar A vector of strings specifying the names of covariates derived from traits. See \code{\link{pheno2covar}}.
    p_covar = NULL,
    #' @field g_covar A vector of strings specifying the names of covariates derived from genetic markers.
    #' See \code{\link{marker2covar}}.
    g_covar = NULL,
    #' @field p_covar_table A matrix holding the individual values for each trait-derived covariate. See \code{\link{pheno2covar}}.
    p_covar_table = NULL,
    #' @field g_covar_table A matrix holding the individual values for each marker-derived covariate. See \code{\link{marker2covar}}.
    g_covar_table = NULL,
    #' @field model_family Indicates the model family of the phenotypes. This can be either "gaussian" or "binomial".
    #' If this argument is length 1, all phenotypes will be assigned to the same family. Phenotypes can be assigned
    #' different model families by providing a vector of the same length as the number of phenotypes, indicating how
    #' each phenotype should be modeled. See \code{\link{singlescan}}.
    model_family = NULL,
    #' @field scan_what A string indicating whether "eigentraits", "normalized.traits", or "raw.traits" should be analyzed.
    #' See \code{\link{get.pheno}}.
    scan_what = NULL,
    #' @field ET A matrix holding the eigentraits to be analyzed.
    ET = NULL,
    #' @field singular_values Added by \code{\link{get.eigentraits}}. A vector holding the singular values from the singular
    #' value decomposition of the trait matrix. They are used in rotating the final direct influences back to trait space
    #' from eigentrait space. See \code{\link{get.eigentraits}} and \code{\link{direct.influence}}.
    singular_values = NULL,
    #' @field right_singular_vectors Added by \code{\link{get.eigentraits}}. A matrix containing the right singular vectors
    #' from the singular value decomposition of the trait matrix. They are used in rotating the final direct influences
    #' back to trait space from eigentrait space. See \code{\link{get.eigentraits}} and \code{\link{direct.influence}}.
    right_singular_vectors = NULL,
    #' @field traits_scaled Whether the traits should be mean-centered and standardized before analyzing.
    traits_scaled = NULL,
    #' @field traits_normalized Whether the traits should be rank Z normalized before analyzing.
    traits_normalized = NULL,
    #' @field var_to_var_influences_perm added in \code{\link{error.prop}}. The list of results from the error propagation
    #' of permuted coefficients.
    var_to_var_influences_perm = NULL,
    #' @field var_to_var_influences added in \code{\link{error.prop}}. The list of results from the error propagation of coefficients.
    var_to_var_influences = NULL,
    #' @field pval_correction Options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none".
    pval_correction = NULL,
    #' @field var_to_var_p_val The final table of cape interaction results calculated by \code{\link{error.prop}}.
    var_to_var_p_val = NULL,
    #' @field max_var_to_pheno_influence The final table of cape direct influences of markers to traits calculated
    #' by \code{\link{direct.influence}}.
    max_var_to_pheno_influence = NULL,
    #' @field full_net An adjacency matrix holding significant cape interactions between individual markers. See
    #' \code{\link{plotNetwork}} and \code{\link{get.network}}.
    full_net = NULL,
    #' @field use_kinship Whether to use a kinship correction in the analysis.
    use_kinship = NULL,
    #' @field kinship_type which type of kinship matrix to use
    kinship_type = NULL,
    #' @field transform_to_phenospace whether to transform to phenospace or not.
    transform_to_phenospace = NULL,
    #' @description
    #' Assigns variables from the parameter file to attributes in the Cape object.
    assign_parameters = function() {
      parameter.table <- read.parameters(self$parameter_file, self$yaml_parameters)
      for(name in names(parameter.table)){
        val <- parameter.table[[name]]
        self[[name]] <- val
      }
    },
    #' @description
    #' Checks the dimensionality of inputs and its consistency.
    check_inputs = function() {
      stopifnot(length(self$chromosome) == length(self$marker_location))
      stopifnot(length(self$chromosome) == length(self$marker_num))
      stopifnot(length(self$chromosome) == length(self$geno_names$locus))
    },
    #' @description
    #' Checks genotype names.
    check_geno_names = function() {
      # TODO make sure that individual names match between the pheno object, geno object, and geno names
      stopifnot(TRUE)
    },
    #' @description
    #' Initialization method.
    #' @param parameter_file string, full path to YAML file with initialization
    #'   parameters
    #' @param yaml_parameters string representing YAML CAPE parameters. See the
    #'   vignette for more descriptions of individual parameter settings.
    #' @param results_path string, full path to directory for storing results
    #'   (optional, a directory will be created if one is not specified)
    #' @param save_results Whether to save cape results. Defaults to TRUE.
    #' @param use_saved_results Whether to use existing results from a 
    #'   previous run. This can save time if re-running an analysis, but
    #'   can lead to problems if the old run and new run have competing settings.
    #'   If errors arise, and use_saved_results is set to TRUE, try setting it 
    #'   to FALSE, or deleting previous results.
    #' @param pheno A matrix containing the traits to be analyzed. Traits are in
    #'    columns and individuals are in rows.
    #' @param chromosome A vector the same length as the number of markers indicating
    #'    which chromosome each marker lives on.
    #' @param marker_num A vector the same length as the number of markers indicating
    #'    the index of each marker
    #' @param marker_location A vector the same length as the number of markers indicating
    #'    the genomic position of each marker. The positions are primarily used for plotting
    #'    and can be in base pairs, centiMorgans, or dummy variables.
    #' @param geno_names The dimnames of the genotype array. The genotype array is a three-dimensional
    #'   array in which rows are individuals, columns are alleles, and the third dimension houses
    #'   the markers. Genotypes are pulled for analysis using \code{\link{get.geno}} based on
    #'   geno_names. Only the individuals, alleles, and markers listed in geno_names are
    #'   taken from the genotype matrix. Functions that remove markers and individuals from
    #'   analysis always operate on geno_names in addition to other relevant slots.
    #'   The names of geno_names must be "mouse", "allele", "locus."
    #' @param geno A three dimentional array holding genotypes for each animal for each allele
    #'   at each marker. The genotypes are continuously valued probabilities ranging from 0 to 1. 
    #'   The dimnames of geno must be "mouse", "allele", and "locus," even if the individuals are
    #'   not mice.
    #' @param .geno_for_pairscan A two-dimensional matrix holding the genotypes that will be analyzed
    #'   in the pairscan. Alleles are in columns and individuals are in rows. As in the geno array, 
    #'   values are continuous probabilities ranging from 0 to 1.
    #' @param peak_density The density parameter for \code{\link{select.markers.for.pairscan}}.
    #'   Determines how densely markers under an individual effect size peak are selected 
    #'   for the pairscan if marker_selection_method is TRUE. Defaults to 0.5.
    #' @param window_size The window size used by \code{\link{select.markers.for.pairscan}}.
    #'   It specifies how many markers are used to smooth effect size curves for automatic peak
    #'   identification. If set to NULL, window_size is determined automatically. Used when 
    #'   marker_selection_method is TRUE.
    #' @param tolerance The wiggle room afforded to \code{\link{select.markers.for.pairscan}} in 
    #'   finding a target number of markers. If num_alleles_in_pairscan is 100 and the tolerance 
    #'   is 5, the algorithm will stop when it identifies anywhere between 95 and 105 markers 
    #'   for the pairscan.
    #' @param ref_allele A string of length 1 indicating which allele to use as the reference allele.
    #'   In two-parent crosses, this is usually allele A. In DO/CC populations, we recommend using
    #'   B as the reference allele. B is the allele from the C57Bl6/J mouse, which is often used as
    #'   a reference strain.
    #' @param alpha The significance level for calculating effect size thresholds in the 
    #'   \code{\link{singlescan}}. If singlescan_perm is 0, this parameter is ignored.
    #' @param covar_table A matrix of covariates with covariates in columns and individuals
    #'   in rows. Must be numeric.
    #' @param num_alleles_in_pairscan The number of alleles to test in the pairwise scan. 
    #'   Because Cape is computationally intensive, we usually need to test only a subset
    #'   of available markers in the pairscan, particularly if the kinship correction is
    #'   being used.
    #' @param max_pair_cor the maximum Pearson correlation between two markers. If their
    #'   correlation exceeds this value, they will not be tested against each other in the
    #'   pairscan. This threshold is set to prevent false positive arising from testing
    #'   highly correlated markers. If this value is set to NULL, min_per_genotype must
    #'   be specified.
    #' @param min_per_genotype minimum The minimum number of individuals allowable per
    #'   genotype combination in the pair scan. If for a given marker pair, one of the 
    #'   genotype combinations is underrepresented, the marker pair is not tested. If 
    #'   this value is NULL, max_pair_cor must be specified.
    #' @param pairscan_null_size The total size of the null distribution.
    #'   This is DIFFERENT than the number of permutations to run. Each permutation
    #'   generates n choose 2 elements for the pairscan. So for example, a permutation
    #'   that tests 100 pairs of markers will generate a null distribution of size 4950.
    #'   This process is repeated until the total null size is reached. If the null size
    #'   is set to 5000, two permutations of 100 markers would be done to get to a null
    #'   distribution size of 5000.
    #' @param p_covar A vector of strings specifying the names of covariates derived
    #'   from traits. See \code{\link{pheno2covar}}.
    #' @param g_covar A vector of strings specifying the names of covariates derived 
    #'   from genetic markers. See \code{\link{marker2covar}}.
    #' @param p_covar_table A matrix holding the individual values for each
    #'   trait-derived covariate. See \code{\link{pheno2covar}}.
    #' @param g_covar_table A matrix holding the individual values for each 
    #'   marker-derived covariate. See \code{\link{marker2covar}}.
    #' @param model_family Indicates the model family of the phenotypes
    #'   This can be either "gaussian" or "binomial". If this argument
    #'   is length 1, all phenotypes will be assigned to the same
    #'   family. Phenotypes can be assigned different model families by
    #'   providing a vector of the same length as the number of phenotypes,
    #'   indicating how each phenotype should be modeled. See \code{\link{singlescan}}.
    #' @param scan_what A string indicating whether "eigentraits", "normalized.traits", or 
    #'   "raw.traits" should be analyzed. See \code{\link{get.pheno}}.
    #' @param ET A matrix holding the eigentraits to be analyzed.
    #' @param singular_values Added by \code{\link{get.eigentraits}}. A vector holding 
    #'   the singular values from the singular
    #'   value decomposition of the trait matrix. They are used in rotating the 
    #'   final direct influences back to trait space from eigentrait space. See
    #'   \code{\link{get.eigentraits}} and \code{\link{direct.influence}}.
    #' @param right_singular_vectors Added by \code{\link{get.eigentraits}}. A matrix 
    #'   containing the right singular vectors from the singular
    #'   value decomposition of the trait matrix. They are used in rotating the 
    #'   final direct influences back to trait space from eigentrait space. See
    #'   \code{\link{get.eigentraits}} and \code{\link{direct.influence}}.
    #' @param traits_scaled Whether the traits should be mean-centered and standardized
    #'   before analyzing.
    #' @param traits_normalized Whether the traits should be rank Z normalized before
    #'   analyzing.
    #' @param var_to_var_influences_perm added in \code{\link{error.prop}} 
    #'  The list of results from the error propagation of permuted coefficients.
    #' @param var_to_var_influences added in \code{\link{error.prop}} 
    #'  The list of results from the error propagation of coefficients.
    #' @param pval_correction Options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
    #' @param var_to_var_p_val The final table of cape interaction results calculated by \code{\link{error.prop}}.
    #' @param max_var_to_pheno_influence The final table of cape direct influences of markers to traits
    #'  calculated by \code{\link{direct.influence}}.
    #' @param full_net An adjacency matrix holding significant cape interactions between
    #'  individual markers. See \code{\link{plotNetwork}} and \code{\link{get.network}}.
    #' @param use_kinship Whether to use a kinship correction in the analysis.
    #' @param kinship_type Which type of kinship matrix to use. Either "overall" or "ltco."
    #' @param transform_to_phenospace whether to transform to phenospace or not.
    initialize = function(
      parameter_file = NULL,
      yaml_parameters = NULL,
      results_path = NULL,
      save_results = TRUE,
      use_saved_results = TRUE,
      pheno = NULL,
      chromosome = NULL,
      marker_num = NULL,
      marker_location = NULL,
      geno_names = NULL,
      geno = NULL,
      .geno_for_pairscan = NULL,
      peak_density = NULL,
      window_size = NULL,
      tolerance = NULL, 
      ref_allele = NULL,
      alpha = NULL,
      covar_table = NULL,
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
      kinship_type = NULL,
      transform_to_phenospace = NULL
    ) {
      self$parameter_file <- parameter_file
      self$yaml_parameters <- yaml_parameters
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
      self$pheno <- pheno
      self$chromosome <- chromosome
      self$marker_num <- marker_num
      self$marker_location <- marker_location
      self$geno <- geno
      self$peak_density <- peak_density
      self$window_size <- window_size
      self$tolerance <- tolerance
      self$geno_names <- geno_names
      if (!missing(ref_allele)) {
        stopifnot(is.character(ref_allele))
        self$ref_allele <- ref_allele
      }
      if (is.null(alpha)) {
        self$alpha <- c(0.05, 0.01)
      } else {
        self$alpha <- alpha 
      }
      self$covar_table <- covar_table
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
      self$kinship_type <- kinship_type
      self$transform_to_phenospace <- transform_to_phenospace
      # assign parameters from the parameter_file
      self$assign_parameters()
      self$check_inputs()
      self$check_geno_names()
      #check.bad.markers(self)
      # TODO make sure that individual names match between the pheno object, geno object, and geno names
    },
    #' @description
    #' Plot Eignentraits
    #' @param filename filename of result plot
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
    #' @description
    #' Plot results of single-locus scans
    #' @param filename filename of result plot.
    #' @param singlescan.obj a singlescan object from \code{\link{singlescan}}
    #' @param width width of result plot, default is 20.
    #' @param height height of result plot, default is 6.
    #' @param units units of result plot, default is "in".
    #' @param res resolution of result plot, default is 300.
    #' @param standardized If TRUE t statistics are plotted. If FALSE, effect sizes are plotted, default is TRUE
    #' @param allele.labels A vector of labels for the alleles if different that those
    #' stored in the data.object.
    #' @param alpha the alpha significance level. Lines for significance values will only
    #' be plotted if n.perm > 0 when \code{\link{singlescan}} was run. And only alpha values
    #' specified in \code{\link{singlescan}} can be plotted.
    #' @param include.covars Whether to include covariates in the plot.
    #' @param line.type as defined in plot
    #' @param pch see the "points()" R function. Default is 16 (a point).
    #' @param cex see the "points()" R function. Default is 0.5.
    #' @param lwd line width, default is 3.
    #' @param traits a vector of trait names to plot. Defaults to all traits.
    plot_singlescan = function(filename, singlescan.obj, width = 20, height = 6, units = "in", res = 300, 
                               standardized = TRUE, allele.labels = NULL, alpha = alpha, include.covars = TRUE, 
                               line.type = "l", pch = 16, cex = 0.5, lwd = 3, traits = NULL) {
      
      full.path <- file.path(self$results_path, filename)
      
      jpeg(full.path, width = width, height = height, units = units, res = res)
      plotSinglescan(self, singlescan.obj = singlescan.obj, standardized = standardized, allele.labels = allele.labels, 
                     alpha = alpha, include.covars = include.covars, line.type = line.type, pch = pch, cex = cex, 
                     lwd = lwd, traits = traits)
      dev.off()
      
    },
    #' @description
    #' Plot the result of the pairwise scan
    #' @param filename filename of result plot.
    #' @param pairscan.obj a pairscan object from \code{\link{pairscan}}
    #' @param phenotype The names of the phenotypes to be plotted. If NULL, all phenotypes are plotted.
    #' @param show.marker.labels If TRUE marker labels are plotted along the axes. If FALSE, they are omitted.
    #' @param show.alleles If TRUE, the allele of each marker is indicated by color.
    plot_pairscan = function(filename, pairscan.obj, phenotype = NULL, 
                             show.marker.labels = TRUE, show.alleles = FALSE) {
      
      # filename is usually "Pairscan.Regression.pdf"
      
      full.path <- file.path(self$results_path, filename)
      
      plotPairscan(self, pairscan.obj, phenotype = phenotype, pdf.label = full.path, 
                   show.marker.labels = show.marker.labels, show.alleles = show.alleles)
      
    },
    #' @description
    #' Plot cape coefficients
    #' @param filename filename of result plot.
    #' @param width width of result plot, default is 10.
    #' @param height height of result plot, default is 7.
    #' @param p.or.q A threshold indicating the maximum p value (or q value if FDR was used) of significant 
    #' interactions and main effects.
    #' @param standardize Whether to plot effect sizes (FALSE) or standardized effect sizes (TRUE), 
    #' default is TRUE.
    #' @param not.tested.col The color to use for marker pairs not tested. Takes the same values as 
    #' pos.col and neg.col, default is "lightgray".
    #' @param covar.width See pheno.width. This is the same effect for covariates.
    #' @param pheno.width Each marker and trait gets one column in the matrix. If there are many markers, 
    #' this makes the effects on the traits difficult  to see. pheno.width increases the number of columns
    #' given to the phenotypes. For example, if pheno.width = 11, the phenotypes will be shown 11 times wider
    #' than individual markers.
    plot_variant_influences = function(filename, width = 10, height = 7,
                                       p.or.q = p.or.q, standardize = FALSE, 
                                       not.tested.col = "lightgray", 
                                       covar.width = NULL, pheno.width = NULL) {
      
      full.path <- file.path(self$results_path, filename)
      
      if (endsWith(full.path, "pdf")) {
        pdf(full.path, width = width, height = height)
      } else if (endsWith(full.path, "jpg")) {
        jpeg(full.path, quality = 100)
      }
      
      plotVariantInfluences(self, p.or.q = p.or.q, standardize = FALSE, 
                            not.tested.col = "lightgray", 
                            covar.width = NULL, pheno.width = NULL)
      dev.off()
      
    },
    #' @description
    #' Plots cape results as a circular network
    #' @param filename filename of result plot.
    #' @param label.gap A numeric value indicating the size of the gap the chromosomes and their labels, 
    #' default is 10.
    #' @param label.cex A numeric value indicating the size of the labels, default is 1.5.
    #' @param show.alleles TRUE show the alleles, FALSE does not show alleles. Default is FALSE.
    plot_network = function(filename, label.gap = 10, label.cex = 1.5, show.alleles = FALSE) {
      
      full.path <- file.path(self$results_path, filename)
      if (endsWith(full.path, "pdf")) {
        pdf(full.path)
      } else if (endsWith(full.path, "jpg")) {
        jpeg(full.path)
      }

      plotNetwork(self, label.gap = label.gap, label.cex = label.cex, show.alleles = show.alleles)
      dev.off()
    },
    #' @description
    #' Plot the final epistatic network in a traditional network view.
    #' @param filename filename of result plot.
    #' @param zoom Allows the user to zoom in and out on the image if the network is either 
    #' running off the edges of the plot or too small in the middle of the plot, default is 1.2.
    #' @param node.radius The size of the pie chart for each node, default is 0.3.
    #' @param label.nodes A logical value indicating whether the nodes should be labeled.
    #' Users may want to remove labels for large networks, default is TRUE.
    #' @param label.offset The amount by which to offset the node labels from the center of
    #' the nodes, default is 0.4.
    #' @param label.cex The size of the node labels, default is 0.5.
    #' @param bg.col The color to be used in pie charts for non-significant main effects.
    #' Takes the same values as pos.col, default is "lightgray".
    #' @param arrow.length The length of the head of the arrow, default is 0.1.
    #' @param layout.matrix Users have the option of providing their own layout matrix for the
    #' network. This should be a two column matrix indicating the x and y coordinates of each 
    #' node in the network, default is "layout_with_kk".
    #' @param legend.position The position of the legend on the plot, default is "topright".
    #' @param edge.lwd The thickness of the arrows showing the interactions, default is 1.
    #' @param legend.radius The size of the legend indicating which pie piece corresponds to which
    #' traits, default is 2.
    #' @param legend.cex The size of the labels in the legend, default is 0.7.
    #' @param xshift A constant by which to shift the x values of all nodes in the network,
    #' default is -1.
    plot_full_network = function(filename, zoom = 1.2, node.radius = 0.3, label.nodes = TRUE, label.offset = 0.4, label.cex = 0.5, 
                                 bg.col = "lightgray", arrow.length = 0.1, layout.matrix = "layout_with_kk", legend.position = "topright", 
                                 edge.lwd = 1, legend.radius = 2, legend.cex = 0.7, xshift = -1) {
      
      full.path <- file.path(self$results_path, filename)
      
      if (endsWith(full.path, "pdf")) {
        pdf(full.path)
      } else if (endsWith(full.path, "jpg")) {
        jpeg(full.path)
      }

      plotFullNetwork(self, zoom = zoom, node.radius = node.radius, label.nodes = label.nodes, label.offset = label.offset, label.cex = label.cex, 
                      bg.col = bg.col, arrow.length = arrow.length, layout.matrix = layout.matrix, legend.position = legend.position, 
                      edge.lwd = edge.lwd, legend.radius = legend.radius, legend.cex = legend.cex, xshift = xshift)
      dev.off()
    },
    #' @description
    #' Write significant cape interactions to a csv file.
    #' @param filename filename of csv file
    #' @param p.or.q A threshold indicating the maximum adjusted p value considered 
    #' significant. If an FDR method has been used to correct for multiple testing, 
    #' this value specifies the maximum q value considered significant, default is 0.05.
    #' @param include.main.effects Whether to include main effects (TRUE) or only
    #' interaction effects (FALSE) in the output table, default is TRUE.
    write_variant_influences = function(filename, p.or.q = 0.05, 
    include.main.effects = TRUE) {
      
      full.path <- file.path(self$results_path, filename)
      
      writeVariantInfluences(self, p.or.q = max(c(p.or.q, 0.2)), 
      	include.main.effects = include.main.effects, filename = full.path)
    },
    #' @description
    #' Set phenotype
    #' @param val phenotype value.
    set_pheno = function(val) {
      self$pheno <- val
      invisible(self)
    },
    #' @description
    #' Set genotype
    #' @param val genotype value.
    set_geno = function(val) {
      self$geno <- val
      invisible(self)
    },
    #' @description
    #' Create covariate table
    #' @param value covariate values
    create_covar_table = function(value) {
      
      marker.locale <- get.col.num(self$pheno, value)
      
      # make a separate covariate table, then modify the dimnames
      # in the genotype object to include the covariates
      # do not modify the genotype object
      
      self$covar_table <- self$pheno[,marker.locale,drop=FALSE]
      rownames(self$covar_table) <- rownames(self$pheno)
      
      # take the phenotypes made into markers out of the phenotype matrix
      self$pheno <- self$pheno[,-marker.locale]
      self$geno_names[[3]] <- c(self$geno_names[[3]], value)
      self$chromosome <- c(self$chromosome, rep(0, length(value)))
      self$marker_location <- c(self$marker_location, 1:length(value))
      invisible(self)
    },
    #' @description
    #' Save to RDS file
    #' @param object data to be saved.
    #' @param filename filename of result RData file.
    save_rds = function(object, filename) {
      # only save the results RData file if save_results is TRUE
      if (self$save_results) {
        full.path <- file.path(self$results_path, filename)
        saveRDS(object, full.path)
      }
    },
    #' @description
    #' Read RDS file
    #' @param filename RData filename to be read.
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
