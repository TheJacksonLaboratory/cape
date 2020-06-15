#' Generate a null distribution for the pairscan.
#' 
#' This script generates a null distribution
#' for the pairscan. For each permutation,
#' it runs a single scan and selects the top
#' N markers. It then uses these markers to
#' perform a permutation of the pairscan.
#' the null distribution generated here uses
#' a fixed number of the TOP ranking markers
#' from the permuted single scan
#' we need to update it to use select.markers.for.pairscan
#' if marker.selection.method is netwas, you need to provide
#' a list of genes from the netWAS analysis
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#' @param scan.what A character string uniquely identifying whether eigentraits
#'   or raw traits should be scanned. Options are "eigentraits", "raw.traits"
#' @param pairscan.null.size
#' @param max.pair.cor A numeric value between 0 and 1 indicating the maximum
#'   Pearson correlation that two markers are allowed. If the correlation
#'   between a pair of markers exceeds this threshold, the pair is not tested.
#'   If this value is set to NULL, min.per.genotype must have a numeric value.
#' @param min.per.geno The minimum number of individuals allowable per
#'   genotype. If for a given marker pair, one of the genotypes is
#'   underrepresented, the marker pair is not tested. If this value is NULL,
#'   max.pair.cor must have a numeric value.
#' @param model.family Indicates the model family of the phenotypes. This can be 
#'   either "gaussian" or "binomial".
#' @param marker.selection.method options are "top.effects", "uniform", "effects.dist", "by.gene"
#' @param run.parallel
#' @param n.cores
#' @param gene.list boolean, only required for "by.gene" marker selection method
#' @param verbose boolean, default = FALSE
#' 
pairscan.null <- function(data.obj, geno.obj = NULL, scan.what = c("eigentraits", "raw.traits"), 
  pairscan.null.size = NULL, max.pair.cor = NULL, min.per.geno = NULL, 
  model.family = "gaussian", 
  marker.selection.method = c("top.effects", "uniform", "effects.dist", "by.gene"), 
  run.parallel = FALSE, gene.list = NULL, n.cores = 4, verbose = FALSE){
  
  marker.selection.method <- data.obj$marker_selection_method
  ref.allele <- data.obj$ref_allele
  
  
  if(is.null(pairscan.null.size)){
    stop("The total number of permutations must be specified.")
  }
  
  
  #If the user does not specify a scan.what, 
  #default to eigentraits, basically, if eigen,
  #et, or ET are anywhere in the string, use the
  #eigentraits, otherwise, use raw phenotypes
  type.choice <- c(grep("eig", scan.what, ignore.case = TRUE), grep("ET", scan.what, ignore.case = TRUE)) #look for any version of eigen or eigentrait, the user might use.
  if(length(type.choice) > 0){ #if we find any, use the eigentrait matrix
    pheno <- data.obj$ET
  }else{
    pheno <- data.obj$pheno #otherwise, use the raw phenotype matrix
  }
  
  num.pheno <- dim(pheno)[2]
  #use the full genotype matrix to select 
  #markers for generating the null in the 
  #pairscan
  
  geno <- get.geno(data.obj, geno.obj)
  
  
  #make a list to hold the results. 
  #Tables from each of the phenotypes will be
  #put into this list
  results.perm.list <- vector(mode = "list", length = num.pheno)
  
  #generate the null distribution for the pairscan 
  #do a singlescan on each permuted trait.
  #find the top n markers
  #combine these into a unique list
  #do the pairscan on these markers for all traits
  
  if(verbose){
    cat("\nGenerating null distribution...\n")
  }
  
  all.pairs.tested <- NULL
  
  n.top.markers <- ncol(data.obj$geno_for_pairscan)
  final.perm <- 1
  while(final.perm < pairscan.null.size){
    perm.order <- sample(1:dim(pheno)[1])
    
    
    if(marker.selection.method != "by.gene" && marker.selection.method != "from.list"){
      single.scan.result <- array(NA, dim = c(length(data.obj$geno_names[[3]]), num.pheno, (dim(geno)[[2]]-1)))
      dimnames(single.scan.result) <- list(data.obj$geno_names[[3]], colnames(pheno), dimnames(geno)[[2]][-which(
        dimnames(geno)[[2]] == ref.allele)])
      
      if(verbose){cat("Performing single marker scans of permuted traits.\n")}
      
      for(p in 1:num.pheno){ 
        if(verbose){cat("\t", colnames(pheno)[p], "...\n")}
        one.singlescan.tstats <- one.singlescanDO(phenotype.vector = pheno[perm.order,p], 
        		genotype.mat = geno, ref.allele = ref.allele, model.family = model.family, 
        		run.parallel = run.parallel, n.cores = n.cores)
        single.scan.result[,p,] <- one.singlescan.tstats
      }
      
      if(verbose){cat("Selecting markers for permuted pairscan...\n")}				
      #use this singlescan to select markers for a permuted pairscan
      
      if(marker.selection.method == "top.effects"){
        perm.data.obj <- select.markers.for.pairscan(data.obj, singlescan.obj = single.scan.result, geno.obj, 
                                                     num.alleles = n.top.markers, peak.density = data.obj$peak_density, 
                                                     window.size = data.obj$window_size, tolerance = data.obj$tolerance, 
                                                     plot.peaks = FALSE, verbose = TRUE)
      }
      if(marker.selection.method == "uniform"){
        perm.data.obj <- select.markers.for.pairscan.uniform(data.obj, geno.obj, num.alleles = ncol(data.obj$geno_for_pairscan), verbose = FALSE)	
      }
      if(marker.selection.method == "effects.dist"){
        perm.data.obj <- select.markers.for.pairscan.dist(data.obj, singlescan.obj = single.scan.result, geno.obj, verbose = FALSE)		
      }
    }else{ 
      
      if(marker.selection.method == "by.gene"){
        #if we are using a gene-based method
        #use a permuted gene list to select
        #SNPs near genes
        perm.data.obj <- select.markers.for.pairscan.by.gene(data.obj, ref.allele = ref.allele, geno.obj = geno.obj, 
                                                             gene.list = sample(gene.list), num.snps = ncol(data.obj$geno_for_pairscan), 
                                                             organism = data.obj$organism)
      }
      if(marker.selection.method == "from.list"){
        single.scan.result <- list("ref.allele" = ref.allele)
        specific.markers <- colnames(data.obj$geno_for_pairscan)
        perm.data.obj <- select.markers.for.pairscan(data.obj, singlescan.obj = single.scan.result, geno.obj, specific.markers = specific.markers)
      }
    }
    
    
    if(verbose){cat("\tGetting markers for permuted pairscan...\n")}
    top.marker.pairs <- get.pairs.for.pairscan(gene = perm.data.obj$geno_for_pairscan, max.pair.cor = max.pair.cor, min.per.genotype = min.per.geno, verbose = FALSE)
    total.pairs <- nrow(top.marker.pairs)
    num.to.add <- 10
    #we don't need to do extra permutations
    #so trim the final pair matrix down to get only
    #the specified number of permutations plus a few
    #because some pairs are always rejected
    
    if(final.perm+dim(top.marker.pairs)[1] > pairscan.null.size){
      num.needed <- pairscan.null.size - final.perm
      #testing just one pair was messing this up, so 
      #always test at least two pairs
      top.marker.pairs <- top.marker.pairs[1:(num.needed+(min(c(num.to.add, total.pairs)))),,drop=FALSE]
    }
    
    if(verbose){cat("\tTesting", dim(top.marker.pairs)[1], "pairs...\n")}
    all.pairs.tested <- rbind(all.pairs.tested, top.marker.pairs)
    
    #run the pairscan for each permuted phenotype and the pairs we just found
    if(verbose){cat("Performing marker pair scans of permuted traits...\n")}
    for(p in 1:num.pheno){
      if(verbose){cat("\t", colnames(pheno)[p], "...\n")}
      #run a pairscan on these markers and each permuted phenotype
      pairscan.results <- one.pairscan.parallel(perm.data.obj, phenotype.vector = pheno[perm.order,p], 
                                                genotype.matrix = perm.data.obj$geno_for_pairscan, 
                                                paired.markers = top.marker.pairs, n.perm = 0, 
                                                run.parallel = run.parallel, n.cores = n.cores, 
                                                verbose = verbose)
      
      #integrate the results into the permutation object
      one.perm <- pairscan.results[[1]]
      #because there will be different numbers of markers each time, just take 
      #the marker names, the intercept, and the effects for marker1 marker2 and 
      #their interaction
      # last.col = dim(one.perm[[1]])[2]
      # take.col <- c(1:3, (last.col-2):last.col)
      if(final.perm == 1){ #if this is the first time through, just copy the results into the results.perm.list
        for(i in 1:2){ #get the effects and the se
          # results.perm.list[[p]][[i]] <- one.perm[[i]][,take.col]
          results.perm.list[[p]][[i]] <- one.perm[[i]]
        }
        results.perm.list[[p]][[3]] <- one.perm[[3]]
      }else{
        if(!is.null(one.perm)){
          for(i in 1:length(one.perm)){
            results.perm.list[[p]][[i]] <- rbind(results.perm.list[[p]][[i]], one.perm[[i]])
          }
        }
      }
    }
    final.perm <- dim(results.perm.list[[1]][[1]])[1] #end looping through phenotypes
    if(verbose){cat("\t", final.perm, " null tests: ", round((final.perm/pairscan.null.size)*100), "%...\n", sep = "")} 
  } #end when we have enough permutations
  
  
  names(results.perm.list) <- colnames(pheno)
  results.list <- list("pairscan.perm" = results.perm.list, "pairs.tested.perm" = all.pairs.tested)
  return(results.list)
  
}
