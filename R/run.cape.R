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
#' @param data.obj the S4 class from \code{\link{Cape}}
#' @param geno.obj the genotype object
#' @param parameter.file a full path string to the YAML file containing configuration parameters
#  TODO the config parameters together with the parameters for the data.obj should be sufficient to reproduce a run of CAPE
#' @param results.dir a full path string to an existing empty directory. An error is thrown if the directory is not empty.
#' @param verbose boolean, output goes to stdout
#' @param run.parallel boolean, if TRUE runs certain parts of the code as parallel blocks
#'
#' @return None, output artifacts are saved to the results.dir directory
#'
#' @export
run.cape <- function(data.obj, geno.obj, parameter.file = "cape.parameters.yml", p.or.q = 0.05, path = ".", results.file = "cross.RData",
                     n.cores = 4, run.singlescan = TRUE, run.pairscan = TRUE, error.prop.coef = TRUE, error.prop.perm = TRUE,
                     initialize.only = FALSE, verbose = TRUE, run.parallel = TRUE){
  
  data.obj <- compare.markers(data.obj, geno.obj)
  
  results.base.name <- gsub(".RData", "", results.file)
  
  parameter.table <- yaml::read_yaml(parameter.file)
  
  data.obj$ref_allele <- parameter.table$ref.allele
  
  #===============================================================
  # figure out how to synchronize get.eigentraits and
  # the calculation of the kinship matrix. They both may
  # remove individuals. Maybe trim the kinship object to match
  # the eigentraits
  #===============================================================
  #if we want to use a kinship correction
  if(as.logical(parameter.table$use.kinship)){
    kin.file <- paste0(results.base.name, "_kinship.RData")
    if(file.exists(kin.file)){
      kin.obj <- readRDS(kin.file)
    }else{ #if there isn't a kinship object already, we need to make one
      geno <- get.geno(data.obj, geno.obj)
      missing.vals <- which(is.na(geno))
      
      if(length(missing.vals) > 0){ #if there are missing values in the genotype matrix, 
        #we need to impute the missing values
        imp.data.file <- paste0(results.base.name, "_data_imputed.RData")
        imp.geno.file <- paste0(results.base.name, "_geno_imputed.RData")
        if(!file.exists(imp.geno.file)){ #if the imputation hasn't been done already
          cat("There are missing values in geno.obj. Running impute.missing.geno...\n")
          geno.imp <- impute.missing.geno(data.obj, geno.obj)
          data.obj <- geno.imp$data.obj
          saveRDS(data.obj, imp.data.file)
          geno.obj <- geno.imp$geno.obj
          saveRDS(geno.obj, imp.geno.file)
        }else{ #if the imputation has been done, read in the imputed genotypes
          data.obj <- readRDS(imp.data.file)
          geno.obj <- readRDS(imp.geno.file)
        }
      } #end case for when there are missing values in the genotype object
      kin.obj <- Kinship(data.obj, geno.obj, type = kinship.type, pop = pop, locus = locus)
      saveRDS(kin.obj, kin.file)
    }
  } else {
    kin.obj <- NULL
  }
  
  if(any(!run.singlescan, !run.pairscan, !error.prop.coef, !error.prop.perm)){
    data.obj <- readRDS(results.file)
  }else{
    
    
    if(verbose){cat("Removing unused markers...\n")}
    data.obj <- remove.unused.markers(data.obj, geno.obj)
    combined.data.obj <- delete.underscore(data.obj, geno.obj)
    
    data.obj <- combined.data.obj$data.obj
    geno.obj <- combined.data.obj$geno.obj
    
    if(!is.null(covariates)){
      data.obj <- pheno2covar(data.obj, covariates)
    }
    if(!is.null(marker.covariates)){
      data.obj <- marker2covar(data.obj, geno.obj, markers = marker.covariates)
    }
    
    data.obj <- select.pheno(data.obj, pheno.which = traits)	
    
    
    if(length(grep("e", scan.what, ignore.case = TRUE)) > 0){
      data.obj <- get.eigentraits(data.obj, scale.pheno = as.logical(traits.scaled), normalize.pheno = as.logical(traits.normalized))
      
      pdf("svd.pdf")
      plotSVD(data.obj, orientation = "vertical")
      dev.off()
      
      data.obj <- select.eigentraits(data.obj, traits.which = eig.which)
      
      if(use.kinship){
        #if individuals were deleted from the phenotype matrix, delete these
        #from the kinship object too
        kin.obj <- remove.kin.ind(data.obj, kin.obj)
        saveRDS(kin.obj, kin.file)
      }
    }
    
    saveRDS(data.obj, results.file)
  }
  
  if(initialize.only){
    return(data.obj)
  }
  
  #===============================================================
  # run singlescan
  #===============================================================
  singlescan.results.file <- paste0(results.base.name, ".singlescan.RData")
  
  if(run.singlescan){
    singlescan.obj <- singlescan(data.obj, geno.obj, kin.obj = kin.obj, n.perm = singlescan.perm, ref.allele = ref.allele, alpha = c(0.01, 0.05), scan.what = scan.what, verbose = verbose, run.parallel = run.parallel, n.cores = n.cores, model.family = "gaussian", overwrite.alert = FALSE)
    saveRDS(singlescan.obj, singlescan.results.file)
    
    
    for(ph in 1:ncol(singlescan.obj$singlescan.effects)){
      jpeg(paste0("Singlescan.", colnames(singlescan.obj$singlescan.effects)[ph], ".Standardized.jpg"), width = 20, height = 6, units = "in", res = 300)
      plotSinglescan(data.obj, singlescan.obj = singlescan.obj, standardized = TRUE, allele.labels = NULL, alpha = c(0.05, 0.01), include.covars = TRUE, line.type = "l", pch = 16, cex = 0.5, lwd = 3, traits = colnames(singlescan.obj$singlescan.effects)[ph])
      dev.off()
    }
    
    for(ph in 1:ncol(singlescan.obj$singlescan.effects)){
      jpeg(paste0("Singlescan.", colnames(singlescan.obj$singlescan.effects)[ph], ".Effects.jpg"), width = 20, height = 6, units = "in", res = 300)			
      plotSinglescan(data.obj, singlescan.obj = singlescan.obj, standardized = FALSE, allele.labels = NULL, alpha = c(0.05, 0.01), include.covars = TRUE, line.type = "l", pch = 16, cex = 0.5, lwd = 3, traits = colnames(singlescan.obj$singlescan.effects)[ph])
      dev.off()
    }
    
    
  }else{
    singlescan.obj <- readRDS(singlescan.results.file)
  }
  
  
  #===============================================================
  # run pairscan
  #===============================================================
  pairscan.file <- paste0(results.base.name, ".pairscan.RData")
  if(run.pairscan){
    if(marker.selection.method == "top.effects"){
      data.obj <- select.markers.for.pairscan(data.obj, singlescan.obj, geno.obj, num.alleles = num.alleles.in.pairscan, peak.density = peak.density, verbose = verbose, plot.peaks = FALSE)
    }
    
    if(marker.selection.method == "from.list"){
      specific.markers <- read.table(SNPfile, sep = "\t", stringsAsFactors = FALSE)
      data.obj <- select.markers.for.pairscan(data.obj, singlescan.obj, geno.obj, specific.markers = specific.markers[,1], verbose = verbose, plot.peaks = FALSE)
    }
    
    
    if(marker.selection.method == "uniform"){
      data.obj <- select.markers.for.pairscan.uniform(data.obj, geno.obj, ref.allele = ref.allele, required.markers = NULL, num.alleles = num.alleles.in.pairscan, verbose = verbose)
    }
    
    
    if(marker.selection.method == "by.gene"){
      gene.list.mat <- read.table("gene.list.txt", sep = "\t", stringsAsFactors = FALSE)		
      gene.list <- gene.list.mat[,1]
      data.obj <- select.markers.for.pairscan.by.gene(data.obj, geno.obj, ref.allele = ref.allele, gene.list = gene.list, num.snps = num.alleles.in.pairscan, bp.buffer = bp.buffer, organism = organism)
    } else {
      gene.list <- NULL
    }
    
    saveRDS(data.obj, results.file)
    
    pairscan.obj <- pairscan(data.obj, geno.obj, scan.what = scan.what, pairscan.null.size = pairscan.null.size, min.per.genotype = min.per.geno, max.pair.cor = max.pair.cor, verbose = verbose, num.pairs.limit = Inf, overwrite.alert = FALSE, run.parallel = run.parallel, n.cores = n.cores, gene.list = gene.list, kin.obj = kin.obj)
    saveRDS(pairscan.obj, pairscan.file)
    
    plotPairscan(data.obj, pairscan.obj, phenotype = NULL, pdf.label = "Pairscan.Regression.pdf", show.marker.labels = TRUE, show.alleles = FALSE)
    
    saveRDS(data.obj, results.file)
  } else {
    pairscan.obj <- readRDS(pairscan.file)
  }
  
  
  #===============================================================
  # run reprametrization
  #===============================================================
  
  if(error.prop.coef){
    data.obj <- error.prop(data.obj, pairscan.obj, perm = FALSE, verbose = verbose, n.cores = n.cores, run.parallel = run.parallel)
    saveRDS(data.obj, results.file)
  }
  
  if(error.prop.perm){	
    data.obj <- error.prop(data.obj, pairscan.obj, perm = TRUE, verbose = verbose, n.cores = n.cores, run.parallel = run.parallel)
    saveRDS(data.obj, results.file)
  }
  
  data.obj <- calc.p(data.obj, pval.correction = pval.correction)
  
  if(length(grep("e", scan.what, ignore.case = TRUE)) > 0){
    transform.to.phenospace <- TRUE
  }else{
    transform.to.phenospace <- FALSE	
  }
  
  data.obj <- direct.influence(data.obj, pairscan.obj, transform.to.phenospace = transform.to.phenospace, verbose = TRUE, pval.correction = pval.correction, save.permutations = TRUE, n.cores = n.cores)
  
  saveRDS(data.obj, results.file)
  
  writeVariantInfluences(data.obj, p.or.q = max(c(p.or.q, 0.2)), filename = "Variant.Influences.csv")
  
  pdf("variant.influences.pdf", width = 10, height = 7)
  plotVariantInfluences(data.obj, p.or.q = p.or.q, standardize = FALSE, not.tested.col = "lightgray", covar.width = 30, pheno.width = 30)
  dev.off()
  
  data.obj <- get.network(data.obj, p.or.q = p.or.q, collapse.linked.markers = FALSE)
  data.obj <- get.network(data.obj, p.or.q = p.or.q, threshold.power = 1, collapse.linked.markers = TRUE, plot.linkage.blocks = FALSE)
  saveRDS(data.obj, results.file)
  
  pdf("Network.Circular.pdf")
  plotNetworkDO(data.obj, label.gap = 10, label.cex = 1.5, show.alleles = FALSE)
  dev.off()
  
  if(dim(geno.obj)[2] == 8){
    pdf("Network.Circular.DO.pdf")
    plotNetworkDO(data.obj, label.gap = 10, label.cex = 1.5, show.alleles = TRUE)
    dev.off()		
  }	
  
  pdf("Network.View.pdf")
  net.layout <- plotNetwork2(data.obj, zoom = 1.2, node.radius = 0.3, label.nodes = TRUE, label.offset = 0.4, label.cex = 0.5, bg.col = "lightgray", arrow.length = 0.1, layout.matrix = "layout_with_kk", legend.position = "topright", edge.lwd = 1, legend.radius = 2, legend.cex = 0.7, xshift = -1)
  dev.off()
  
  
  saveRDS(data.obj, results.file)
  
  invisible(data.obj)
  
  
}
