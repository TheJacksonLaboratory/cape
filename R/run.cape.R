#This function runs a cape analysis from a parameter file
#The input is the initial cape data object and a parameter
#file name
#
#and functions loaded 
#if using gene.based marker selection, there must be a file in
#the working directory called gene.list.txt with an ordered list
#of genes in a column
#kinship.type can be either "overall" or "LTCO"
# parameter.file = "cape.parameters.txt"; p.or.q = 0.05; results.file = "cross.RData"; n.cores = 4; run.singlescan = TRUE; run.pairscan = TRUE; error.prop.coef = TRUE; error.prop.perm = TRUE; initialize.only = FALSE; verbose = TRUE; run.parallel = TRUE
#' Runs the CAPE algorithm
#'
#' This function assumes you already have all required libraries and functions loaded.
#'
#' @param pheno.obj the phenotype object
#' @param geno.obj the genotype object
#' @param results.file the name of the saved data.obj RData file. The base name is used as the base name for all saved RData files.
#' @param p.or.q A threshold indicating the maximum adjusted p value considered 
#' @param n.cores integer, default is 4
#' @param run.singlescan boolean, defaul: TRUE
#' @param run.pairscan boolean, default: TRUE
#' @param error.prop.coef, boolean, default: TRUE
#' @param error.prop.perm, boolean, default: TRUE
#' @param initialize.only, boolean, default: FALSE
#' @param verbose boolean, output goes to stdout
#' @param run.parallel boolean, if TRUE runs certain parts of the code as parallel blocks
#' @param param.file yaml parameter file full path
#' @param yaml.params yaml string containing the parameters. Either the param.file or yaml.params can be null.
#' @param results.path paths of the results
#'
#' @return None, output artifacts are saved to the data.obj$results_path directory
#'
#' @export
run.cape <- function(pheno.obj, geno.obj, 
  results.file = "cross.RData", p.or.q = 0.05, n.cores = 4, 
  run.singlescan = TRUE, run.pairscan = TRUE, error.prop.coef = TRUE, 
  error.prop.perm = TRUE, initialize.only = FALSE, verbose = TRUE, 
  run.parallel = FALSE, param.file = NULL, yaml.params = NULL,
  results.path = NULL){
  
  # Instantiate the Cape R6 object
  data.obj <- Cape$new(
  		parameter_file = param.file,
		yaml_parameters = yaml.params,
  		results_path = results.path,
  		pheno = pheno.obj$pheno,
  		chromosome = pheno.obj$chromosome,
  		marker_num = pheno.obj$marker_num,
  		marker_location = pheno.obj$marker_location,
  		geno_names = pheno.obj$geno_names,
  		geno = geno.obj,
  		use_kinship = TRUE
  )
  
  results.base.name <- gsub(".RData", "", results.file)
    
  # since this is the main data.obj, we can't allow it to return FALSE, 
  #check for the file first
  prior.data.obj <- data.obj$read_rds(results.file)
  if (isFALSE(prior.data.obj)) {
    data.obj <- compare.markers(data.obj, geno.obj)
  } else {
    # things can get pretty confusing if these values don't match between 
    #the parameter file and the old data.obj
    prior.data.obj$save_results <- data.obj$save_results
    prior.data.obj$use_saved_results <- data.obj$use_saved_results
    data.obj <- prior.data.obj
  }
  
  #===============================================================
  # figure out how to synchronize get.eigentraits and
  # the calculation of the kinship matrix. They both may
  # remove individuals. Maybe trim the kinship object to match
  # the eigentraits
  #===============================================================
  #if we want to use a kinship correction
  if(as.logical(data.obj$use_kinship)){
    kin.file.name <- paste0(results.base.name, "_kinship.RData")
    kin.obj <- data.obj$read_rds(kin.file.name)
    
    if (isFALSE(kin.obj)) {
      #if there isn't a kinship object already, we need to make one
      kin.obj <- Kinship(data.obj, geno.obj, type = data.obj$kinship_type, 
      pop = data.obj$pop)
      data.obj$save_rds(kin.obj, kin.file.name)
    }

    #===============================================================
    # We need a complete genotype matrix to calculate the kinship
    # adjusted genotypes later on.
    # Check for missing values in the genotype matrix.
    # If there are missing values, impute them.
    # Write out the imputed matrix, or read this in if it already
    # exists.
    #===============================================================
    #we need to impute the missing values
    imp.data.file <- paste0(results.base.name, "_data_imputed.RData")
    imp.geno.file <- paste0(results.base.name, "_geno_imputed.RData")

    # check if there is already a saved genotype object
    geno <- data.obj$read_rds(imp.geno.file)

    if (isFALSE(geno)) {  #if the imputation hasn't been done already
      geno <- get.geno(data.obj, geno.obj)
      missing.vals <- which(is.na(geno))

      if (length(missing.vals) > 0) { #if there are missing values, impute them
        cat("There are missing values in geno.obj. Running impute.missing.geno...\n")
        geno.imp <- impute.missing.geno(data.obj, geno.obj = geno.obj, k = 10, 
        	ind.missing.thresh = 0, marker.missing.thresh = 0, prioritize = "fewer",
        	max.region.size = NULL, min.region.size = NULL, run.parallel = run.parallel,
        	verbose = verbose, n.cores = n.cores)

        # update and save the data.obj
        data.obj <- geno.imp$data.obj
        data.obj$save_rds(data.obj, imp.data.file)

        # update and save the geno.obj
        geno.obj <- geno.imp$geno.obj
        data.obj$save_rds(geno.obj, imp.geno.file)
      
        # recalculate the kinship matrix with the updated objects
        kin.obj <- Kinship(data.obj, geno.obj, type = data.obj$kinship_type, pop = data.obj$pop)
        data.obj$save_rds(kin.obj, kin.file.name)

      } #end case for when there are missing values but no imputed genotypes

    } else { #if the imputation has been done, then it must have been done for the data.obj too
      #data.obj <- geno.imp$data.obj
      geno.obj <- geno
    }

  }
    
  if(verbose){cat("Removing unused markers...\n")}
  data.obj <- remove.unused.markers(data.obj, geno.obj)
  combined.data.obj <- delete_underscore(data.obj, geno.obj)

  data.obj <- combined.data.obj$data.obj
  geno.obj <- combined.data.obj$geno.obj
  
  #because the genotype object can be changed by the above step, 
  #save the final version. (or change the above step so it doesn't 
  final.geno.file <- paste0(results.base.name, "_geno.RData")
  data.obj$save_rds(geno.obj, final.geno.file)

  #str(data.obj$geno_names)
  #str(dimnames(geno.obj))

  if(!is.null(data.obj$covariates)){
    data.obj <- pheno2covar(data.obj, data.obj$covariates)
  }
  if(!is.null(data.obj$marker_covariates)){
    data.obj <- marker2covar(data.obj, geno.obj, markers = data.obj$marker_covariates)
  }
  
  data.obj <- select.pheno(data.obj, pheno.which = data.obj$traits)	
  
  if(length(grep("e", data.obj$scan_what, ignore.case = TRUE)) > 0){
    data.obj <- get.eigentraits(
      data.obj, 
      scale.pheno = as.logical(data.obj$traits_scaled), 
      normalize.pheno = as.logical(data.obj$traits_normalized)
    )
    
    data.obj$plot_svd("svd.pdf")
    data.obj$plot_svd("svd.jpg")
    
    # TODO update select.eigentraits
    data.obj <- select.eigentraits(data.obj, traits.which = data.obj$eig_which)
  }
  
  data.obj$save_rds(data.obj, results.file)

  if(initialize.only){
    return(data.obj)
  }
  
  #===============================================================
  # run singlescan
  #===============================================================
  singlescan.results.file <- paste0(results.base.name, ".singlescan.RData")
  
  singlescan.obj <- data.obj$read_rds(singlescan.results.file)
  
  if (!exists("kin.obj")) {
    kin.obj <- NULL
  }
  
  if (isFALSE(singlescan.obj)) {
    
    if (run.singlescan) {
      
      singlescan.obj <- singlescan(
        data.obj, geno.obj, kin.obj = kin.obj, n.perm = data.obj$singlescan_perm,
        alpha = c(0.01, 0.05), verbose = verbose, run.parallel = run.parallel,
        n.cores = n.cores, model.family = "gaussian", overwrite.alert = FALSE
      )
      
      data.obj$save_rds(singlescan.obj, singlescan.results.file)
      
      for(ph in 1:ncol(singlescan.obj$singlescan.effects)){
        filename <- paste0("Singlescan.", colnames(singlescan.obj$singlescan.effects)[ph], ".Standardized.jpg")
        data.obj$plot_singlescan(filename, singlescan.obj, width = 20, height = 6, 
          units = "in", res = 300, standardized = TRUE, allele.labels = NULL, 
          alpha = c(0.05, 0.01), include.covars = TRUE, line.type = "l", pch = 16, cex = 0.5, 
          lwd = 3, traits = colnames(singlescan.obj$singlescan.effects)[ph])
      }
      
      for(ph in 1:ncol(singlescan.obj$singlescan.effects)){
        filename <- paste0("Singlescan.", colnames(singlescan.obj$singlescan.effects)[ph], ".Effects.jpg")
        data.obj$plot_singlescan(filename, singlescan.obj, width = 20, height = 6, units = "in", res = 300,
          standardized = FALSE, allele.labels = NULL, alpha = c(0.05, 0.01), include.covars = TRUE,
          line.type = "l", pch = 16, cex = 0.5, lwd = 3, traits = colnames(singlescan.obj$singlescan.effects)[ph])
      }
    }
  }
  
  #===============================================================
  # run pairscan
  #===============================================================
  pairscan.file <- paste0(results.base.name, ".pairscan.RData")
  
  pairscan.obj <- data.obj$read_rds(pairscan.file)
  
  if (isFALSE(pairscan.obj) | is.null(data.obj$geno_for_pairscan)) {
    
    if (run.pairscan) {
      
      marker.selection.method <- data.obj$marker_selection_method
      num.alleles.in.pairscan <- data.obj$num_alleles_in_pairscan
      peak.density <- data.obj$peak_density
      max.pair.cor <- data.obj$max_pair_cor
      min.per.genotype <- data.obj$min_per_genotype
      pairscan.null.size <- data.obj$pairscan_null_size
      scan.what <- data.obj$scan_what
      if(marker.selection.method == "top.effects"){
        data.obj <- select.markers.for.pairscan(data.obj, singlescan.obj, geno.obj,
          num.alleles = num.alleles.in.pairscan, peak.density = peak.density,
          verbose = verbose, plot.peaks = FALSE)
      }
      
      if(marker.selection.method == "from.list"){
        snp.file <- file.path(results.path, data.obj$snp_file)
        specific.markers <- as.matrix(read.table(snp.file, sep = "\t", stringsAsFactors = FALSE))
        data.obj <- select.markers.for.pairscan(data.obj, singlescan.obj, geno.obj,
          specific.markers = specific.markers[,1], verbose = verbose, plot.peaks = FALSE)
      }
      
      if(marker.selection.method == "uniform"){
        data.obj <- select.markers.for.pairscan.uniform(data.obj, geno.obj, 
        required.markers = NULL, num.alleles = num.alleles.in.pairscan, verbose = verbose)
      }
      
      if(marker.selection.method == "by.gene"){
        gene.list.mat <- read.table("gene.list.txt", sep = "\t", stringsAsFactors = FALSE)		
        gene.list <- gene.list.mat[,1]
        data.obj <- select.markers.for.pairscan.by.gene(data.obj, geno.obj, gene.list = gene.list, 
                                                        bp.buffer = data.obj$bp_buffer, organism = data.obj$organism)
      } else {
        gene.list <- NULL
      }
      
      data.obj$save_rds(data.obj, results.file)
      
      pairscan.obj <- pairscan(data.obj, geno.obj, scan.what = scan.what, 
        pairscan.null.size = pairscan.null.size, min.per.genotype = min.per.genotype, 
        max.pair.cor = max.pair.cor, verbose = verbose, num.pairs.limit = Inf, 
        overwrite.alert = FALSE, run.parallel = run.parallel, n.cores = n.cores, 
        gene.list = gene.list, kin.obj = kin.obj)
      
      data.obj$save_rds(pairscan.obj, pairscan.file)
      
      data.obj$plot_pairscan("Pairscan.Regression.pdf", pairscan.obj, 
        phenotype = NULL, show.marker.labels = TRUE, show.alleles = FALSE)
      data.obj$plot_pairscan("Pairscan.Regression.jpg", pairscan.obj, 
        phenotype = NULL, show.marker.labels = TRUE, show.alleles = FALSE)
      data.obj$save_rds(data.obj, results.file)
    } 
  }
  
  #===============================================================
  # run reprametrization
  #===============================================================
  
  if(error.prop.coef){
    data.obj <- error.prop(data.obj, pairscan.obj, perm = FALSE, verbose = verbose,
      n.cores = n.cores, run.parallel = run.parallel)
    data.obj$save_rds(data.obj, results.file)
  }
  
  if(error.prop.perm){	
    data.obj <- error.prop(data.obj, pairscan.obj, perm = TRUE, verbose = verbose,
      n.cores = n.cores, run.parallel = run.parallel)
    data.obj$save_rds(data.obj, results.file)
  }
  
  data.obj <- calc.p(data.obj, pval.correction = data.obj$pval_correction)
  
  if(length(grep("e", data.obj$scan_what, ignore.case = TRUE)) > 0){
    transform.to.phenospace <- TRUE
  }else{
    transform.to.phenospace <- FALSE	
  }
  
  data.obj <- direct.influence(data.obj, pairscan.obj, 
    transform.to.phenospace = transform.to.phenospace, verbose = TRUE, 
    pval.correction = data.obj$pval_correction, save.permutations = TRUE, n.cores = n.cores)
  
  data.obj$save_rds(data.obj, results.file)
  
  data.obj$write_variant_influences("Variant.Influences.csv", p.or.q = max(c(p.or.q, 0.2)))

  data.obj$write_variant_influences("Variant.Influences.Interactions.csv", 
    include.main.effects = FALSE, p.or.q = max(c(p.or.q, 0.2)))
  
  data.obj$plot_variant_influences("variant.influences.pdf", width = 10, height = 7,
    p.or.q = p.or.q, standardize = FALSE, not.tested.col = "lightgray", 
    covar.width = NULL, pheno.width = NULL)

  data.obj$plot_variant_influences("variant.influences.jpg", width = 10, height = 7,
    p.or.q = p.or.q, standardize = FALSE, not.tested.col = "lightgray", 
    covar.width = NULL, pheno.width = NULL)

  data.obj <- get.network(data.obj, geno.obj, p.or.q = p.or.q, collapse.linked.markers = FALSE)
  data.obj <- get.network(data.obj, geno.obj, p.or.q = p.or.q, threshold.power = 1, 
    collapse.linked.markers = TRUE, plot.linkage.blocks = FALSE)
  
  data.obj$save_rds(data.obj, results.file)
  
  data.obj$plot_network_do("Network.Circular.pdf", label.gap = 10, label.cex = 1.5, show.alleles = FALSE)
  data.obj$plot_network_do("Network.Circular.jpg", label.gap = 10, label.cex = 1.5, show.alleles = FALSE)

  if(dim(geno.obj)[2] == 8){
    data.obj$plot_network_do("Network.Circular.DO.pdf", label.gap = 10, label.cex = 1.5, show.alleles = TRUE)
    data.obj$plot_network_do("Network.Circular.DO.jpg", label.gap = 10, label.cex = 1.5, show.alleles = TRUE)
  }	
  
  data.obj$plot_full_network("Network.View.pdf", zoom = 1.2, node.radius = 0.3, 
    label.nodes = TRUE, label.offset = 0.4, label.cex = 0.5, bg.col = "lightgray", 
    arrow.length = 0.1, layout.matrix = "layout_with_kk", legend.position = "topright", 
    edge.lwd = 1, legend.radius = 2, legend.cex = 0.7, xshift = -1)
  
  data.obj$plot_full_network("Network.View.jpg", zoom = 1.2, node.radius = 0.3, 
    label.nodes = TRUE, label.offset = 0.4, label.cex = 0.5, bg.col = "lightgray", 
    arrow.length = 0.1, layout.matrix = "layout_with_kk", legend.position = "topright", 
    edge.lwd = 1, legend.radius = 2, legend.cex = 0.7, xshift = -1)
  
  data.obj$save_rds(data.obj, results.file)
  
  invisible(data.obj)
  
  
}
