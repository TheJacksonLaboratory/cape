#' Perform regressions for all pairs of markers and all phenotypes.
#'
#' This function performs the pairwise regression on all selected marker pairs.
#' The phenotypes used can be either eigentraits or raw phenotypes. Permutation
#' testing is also performed.
#'
#' This script performs the pairwise scan on all markers It takes in the data as
#' a cross object. The user has the choice to scan the eigentraits (default) or
#' the original phenotypes. This script also calls the function to do
#' permutations on the 2D scan. It adds the genome-wide threshold for the 2D
#' scan to the data object n.top.markers is used in generating the null. A
#' permutation of the singlescan is run, and the n top markers are used in a
#' permutation of the pairscan. if n.top.markers is NULL, it defaults to the
#' number of markers in geno.for.pairscan
#'
#' @details Not all marker pairs are necessarily tested. Before markers are
#'   tested for interaction, they are checked for several conditions. Pairs are
#'   discarded if (1) at least one of the markers is on the X chromosome, or (2)
#'   there are fewer than min.per.genotype individuals in any of the genotype
#'   combinations.
#'
#'   This function adds an element to the data object reporting the results of
#'   the pair-wise scan: \item{pairscan.results}{The results of the pairwise
#'   scan on the provided phenotype and genotypes.}
#'
#'   If permutations have been performed (n.perm > 0), an additional element is
#'   added to the object reporting the results of the permutation tests:
#'   \item{pairscan.perm}{The results of the permutations of the pairwise scan
#'   on the provided phenotype and genotypes.}
#'
#'   Each of these results elements is itself a list of 3 elements:
#'   \item{pairscan.effects}{A table of effects of each marker pair. The columns
#'   list the effects in the following order: marker1, marker2, the variance of
#'   marker1, the covariance of marker1 and marker2, the variance of marker2,
#'   the covariance of marker1 and the interaction effect, the covariance
#'   between marker2 and ther interaction effect, and the variance of the
#'   interaction.} \item{pairscan.se}{A table of the standard errors from the
#'   test on each marker pair. The columns are identical to those described for
#'   pairscan.effects} \item{model.covariance}{This is a table in which each row
#'   is the linearized matrix of the variance-covariance matrix of each pairwise
#'   regression.}
#'
#'   The results element for the permutation tests has the same structure as for
#'   the pairwise scan except that each row represents the results of one
#'   permutation.
#'
#' @references Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use
#'   of pleiotropy to model genetic interactions in a population. PLoS genetics,
#'   8(10), e1003010. doi:10.1371/journal.pgen.1003010
#'   
#' @seealso \code{\link{select.markers.for.pairscan}}, \code{\link{plotPairscan}}
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#' @param scan.what A character string uniquely identifying whether eigentraits
#'   or raw traits should be scanned. Options are "eigentraits", "raw.traits"
#' @param pairscan.null.size
#' @param max.pair.cor A numeric value between 0 and 1 indicating the maximum
#'   Pearson correlation that two markers are allowed. If the correlation
#'   between a pair of markers exceeds this threshold, the pair is not tested.
#'   If this value is set to NULL, min.per.genotype must have a numeric value.
#' @param min.per.genotype The minimum number of individuals allowable per
#'   genotype. If for a given marker pair, one of the genotypes is
#'   underrepresented, the marker pair is not tested. If this value is NULL,
#'   max.pair.cor must have a numeric value.
#' @param kin.obj a kinship object
#' @param num.pairs.limit A number indicating the maximum number of pairs to
#'   scan. If the number of pairs exceeds this threshold, the function asks for
#'   confirmation before proceeding with the pairwise scan.
#' @param num.perm.limit A number indicating the maximum number of total
#'   permutations that will be performed. If the number of total permutations
#'   exceeds this threshold, the function asks for confirmation before
#'   proceeding with the pairwise scan.
#' @param overwrite.alert
#' @param run.parallel
#' @param n.cores
#' @param gene.list
#' @param verbose boolean, default = FALSE
#'
#' @export
pairscan <- function(data.obj, geno.obj = NULL,
                     scan.what = c("eigentraits", "raw.traits"), pairscan.null.size = NULL, 
                     max.pair.cor = NULL, min.per.genotype = NULL, kin.obj = NULL, 
                     num.pairs.limit = 1e6, num.perm.limit = 1e7, overwrite.alert = TRUE, 
                     run.parallel = TRUE, n.cores = 4, gene.list = NULL, verbose = FALSE) {
  
  marker.selection.method <- data.obj$marker_selection_method
  
  if(!run.parallel){n.cores = 1}
  
  if(is.null(kin.obj)){use.kinship = FALSE}
  if(!is.null(kin.obj)){use.kinship = TRUE}
  data.obj$use_kinship <- use.kinship
  
  
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
  
  cat("Performing pairwise tests...\n")
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
                                         max.pair.cor = max.pair.cor, min.per.genotype, 
                                         verbose = verbose, marker.selection.method = marker.selection.method,
                                         run.parallel = run.parallel, n.cores = n.cores,
                                         gene.list = gene.list)			
    }else{
      pairscan.perm <- pairscan.null(data.obj, geno.obj, scan.what = scan.what, 
                                     pairscan.null.size = pairscan.null.size, 
                                     max.pair.cor = max.pair.cor, min.per.genotype, verbose = verbose, 
                                     marker.selection.method = marker.selection.method, 
                                     run.parallel = run.parallel, n.cores = n.cores,
                                     gene.list = gene.list)
    }
    #add the results to the data object
    pairscan.obj$pairscan.perm <- pairscan.perm$pairscan.perm 
    #add the results to the data object
    pairscan.obj$pairs.tested.perm <- pairscan.perm$pairs.tested.perm 
    
  }
  
  
  return(pairscan.obj) #and return it
  
  }
