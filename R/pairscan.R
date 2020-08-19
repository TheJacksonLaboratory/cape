#' This function performs the pairwise scan on all markers.
#'
#' This function performs the pairwise regression on all selected marker pairs.
#' The phenotypes used can be either eigentraits or raw phenotypes. Permutation
#' testing is also performed.
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#' @param scan.what A character string uniquely identifying whether eigentraits
#'   or raw traits should be scanned. Options are "eigentraits", "raw.traits"
#' @param pairscan.null.size The total size of the null distribution.
#' This is DIFFERENT than the number of permutations to run. Each permutation
#' generates n choose 2 elements for the pairscan. So for example, a permutation
#' that tests 100 pairs of markers will generate a null distribution of size 4950.
#' This process is repeated until the total null size is reached. If the null size
#' is set to 5000, two permutations of 100 markers would be done to get to a null
#' distribution size of 5000.
#' @param max.pair.cor A numeric value between 0 and 1 indicating the maximum
#'   Pearson correlation that two markers are allowed. If the correlation
#'   between a pair of markers exceeds this threshold, the pair is not tested.
#'   If this value is set to NULL, min.per.genotype must have a numeric value.
#' @param min.per.genotype The minimum number of individuals allowable per
#'   genotype combination. If for a given marker pair, one of the genotype combinations is
#'   underrepresented, the marker pair is not tested. If this value is NULL,
#'   max.pair.cor must have a numeric value.
#' @param kin.obj a kinship object calculated by \link{\code{Kinship}}.
#' @param num.pairs.limit A number indicating the maximum number of pairs to
#'   scan. If the number of pairs exceeds this threshold, the function asks for
#'   confirmation before proceeding with the pairwise scan.
#' @param num.perm.limit A number indicating the maximum number of total
#'   permutations that will be performed. If the number of total permutations
#'   exceeds this threshold, the function asks for confirmation before
#'   proceeding with the pairwise scan.
#' @param overwrite.alert If TRUE raises a warning to users not to overwrite 
#'   their data object with a singlescan object. A warning necessary after a 
#'   new version of cape began separating results from different functions into
#'   different results objects
#' @param run.parallel Whether to run the analysis on parallel CPUs
#' @param n.cores The number of CPUs to use if run.parallel is TRUE
#' @param verbose Whether to write progress to the screen
#'
#'
#' @details Not all marker pairs are necessarily tested. Before markers are
#'   tested for interaction, they are checked for several conditions. Pairs are
#'   discarded if (1) at least one of the markers is on the X chromosome, or (2)
#'   there are fewer than min.per.genotype individuals in any of the genotype
#'   combinations.
#'
#' @return This function returns an object assigned to pairscan.obj in 
#' \link{\code{run.cape}}.
#'
#' The results object is a list of five elements:
#' ref.allele: The allele used as the reference for the tests.
#' max.pair.cor: The maximum pairwise correlation between marker pairs
#' pairscan.results: A list with one element per trait. The element for
#' each trait is a list of the following three elements:
#'    pairscan.effects: the effect sizes from the linear models
#'    pairscan.se: the standard erros from the linear models
#'    model.covariance: the model covariance from the linear models.
#' pairscan.perm: The same structure as pairscan.results, but for the
#' permuted data.
#' pairs.tested.perm: A matrix of the marker pairs used in the permutation
#' tests.
#'   
#' @seealso \code{\link{select.markers.for.pairscan}}, \code{\link{plotPairscan}}
#'
#' @export
pairscan <- function(data.obj, geno.obj = NULL,
  scan.what = c("eigentraits", "raw.traits"), pairscan.null.size = NULL, 
  max.pair.cor = NULL, min.per.genotype = NULL, kin.obj = NULL, 
  num.pairs.limit = 1e6, num.perm.limit = 1e7, overwrite.alert = TRUE, 
  run.parallel = FALSE, n.cores = 4, verbose = FALSE) {
  
  marker.selection.method <- data.obj$marker_selection_method
  
  if(!run.parallel){n.cores = 1}
  
  use.kinship <- data.obj$use_kinship
  
  
  if(overwrite.alert){
    choice <- readline(prompt = "Please make sure you assign the output 
                       of this function to a pairscan.obj, and NOT the data.obj. It will 
                       overwrite the data.obj.\nDo you want to continue (y/n) ")
    if(choice == "n"){stop()}
  }
  
  
  pairscan.obj <- list()
  pairscan.obj$ref.allele <- data.obj$ref_allele
  pairscan.obj$max.pair.cor <- data.obj$max_pair_cor
  pairscan.obj$min.per.genotype <- data.obj$min_per_genotype
  
  if(is.null(pairscan.null.size)){
    stop("The final size of the null distribution must be specified.")
  }
  
  pheno <- get.pheno(data.obj, scan.what)	
  
  covar.info <- get.covar(data.obj)
  covar.names <- covar.info$covar.names
  covar.table <- covar.info$covar.table
  
  #find the phenotypic covariates. These will
  #be tested separately, and not as part of a
  #chromosome
  
  if(is.null(data.obj$geno_for_pairscan)){
    stop("select.markers.for.pairscan() must be run before pairscan()")
  }
  
  
  #add the covariates (geno and pheno) 
  #to the genotype matrix so that we 
  #test all pairs
  gene <- get.geno.with.covar(data.obj, geno.obj, g.covar = TRUE, p.covar = TRUE, 
    for.pairscan = TRUE)	
  
  #fill in a matrix to index the marker pairs
  if(verbose){cat("Getting marker pairs for pairscan...\n")}
  pared.marker.mat <- get.pairs.for.pairscan(
    gene,
    covar.names,
    max.pair.cor,
    min.per.genotype,
    run.parallel = run.parallel,
    n.cores = n.cores,
    verbose = verbose
  )
  
  num.pairs <- dim(pared.marker.mat)[1]
  
  if(num.pairs == 0){
    stop("There are no pairs to test. Try raising max.pair.cor or reducing 
         min.per.genotype.")
  }
  
  if(!is.null(num.pairs.limit) && num.pairs > num.pairs.limit){
    cat("\nThe number of pairs (",num.pairs,") exceeds ", num.pairs.limit, ".\n", sep = "")
    go.on <- readline(prompt = "Do you want to continue (y/n)?\n")
    if(length(grep("n", go.on))){
      cat("Stopping pairwise scan...\n")
      return(pairscan.obj)
    }else{
      cat("Continuing pairwise scan...\n")
    }
  }
  
  if(verbose){cat("Performing pairwise tests...\n")}
  #run one.pairscan for each phenotype with results in scanone.result
  if(!use.kinship){
    pairscan.results <- pairscan.noKin(data.obj, pheno.mat = pheno, 
      geno.mat = gene, covar.table = covar.table, marker.pairs = pared.marker.mat, 
      n.perm = pairscan.null.size, verbose = verbose, n.cores = n.cores, 
      run.parallel = run.parallel)
  }else{
    pairscan.results <- pairscan.kin(data.obj, geno.obj = geno.obj, 
      scan.what = scan.what, marker.pairs = pared.marker.mat, kin.obj = kin.obj, 
      verbose = verbose, run.parallel = run.parallel, n.cores = n.cores)
  }	
  
  # print(str(pairscan.results))
  
  pairscan.obj$pairscan.results <- pairscan.results	
  
  if(pairscan.null.size > 0){	
    if(use.kinship){
      pairscan.perm <- pairscan.null.kin(data.obj, geno.obj, kin.obj, 
        scan.what = scan.what, pairscan.null.size = pairscan.null.size, 
        max.pair.cor = max.pair.cor, min.per.genotype, verbose = verbose, 
        marker.selection.method = marker.selection.method, 
        run.parallel = run.parallel, n.cores = n.cores)
    }else{
      pairscan.perm <- pairscan.null(data.obj, geno.obj, scan.what = scan.what, 
        pairscan.null.size = pairscan.null.size, max.pair.cor = max.pair.cor, 
        min.per.genotype, verbose = verbose, marker.selection.method = marker.selection.method, 
        run.parallel = run.parallel, n.cores = n.cores)
    }
    #add the results to the data object
    pairscan.obj$pairscan.perm <- pairscan.perm$pairscan.perm 
    #add the results to the data object
    pairscan.obj$pairs.tested.perm <- pairscan.perm$pairs.tested.perm 
    
  }
  
  
  return(pairscan.obj) #and return it
  
  }
