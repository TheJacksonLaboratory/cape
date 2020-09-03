#' Run the pairscan with a kinship correction
#' 
#' This function is called by \link{\code{pairscan}}
#' when a kinship correction is requested. It adjusts 
#' each variable according to the kinship matrix using
#' \link{\code{kinship.on.the.fly}} and then fits linear
#' pairwise models to the adjusted data.
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#' @param scan.what A character string uniquely identifying whether eigentraits
#'   or raw traits should be scanned. Options are "eigentraits", "raw.traits"
#' @param marker.pairs A two-column matrix containing the marker pairs
#' to be tested.
#' @param kin.obj a kinship object
#' @param verbose A logical value indicating whether to 
#' print progress to the screen
#' @param run.parallel A logical value indicating 
#' whether parallel processing should be used
#' @param n.cores The number of cores to be used if run.parallel is TRUE
#' 
#' @return This function returns a list with three elements. 
#' The elements contain the marker pair effect sizes, the marker
#' pair standard errors, and the covariance matrix for each test.
#' The output is then further processed by \link{\code{pairscan}}.
#' 
pairscan.kin <- function(data.obj, geno.obj, scan.what, marker.pairs, 
kin.obj, verbose = FALSE, run.parallel = FALSE, n.cores = 2){
  
  m = NULL #for appeasing R CMD check
  
  pheno <- get.pheno(data.obj, scan.what)	
  geno <- get.geno.with.covar(data.obj, geno.obj, for.pairscan = TRUE)
  num.pheno <- dim(pheno)[2]
  results.list <- vector(mode = "list", length = num.pheno)
  names(results.list) <- colnames(pheno)
  
  covar.info <- get.covar(data.obj)
  is.char <- as.logical(is.na(suppressWarnings(as.numeric(marker.pairs[1,1]))))
  if(is.char){
    covar.names <- get.marker.name(data.obj, covar.info$covar.names)
  }else{
    covar.names <- get.marker.num(data.obj, covar.info$covar.names)	
  }
  
  p.covar.locale <- which(covar.info$covar.type == "p")
  num.p.covar <- length(p.covar.locale)
  #============================================================================
  #internal functions
  #============================================================================
  
  get.marker.pair.stats <- function(m, kin.dat){
    #if(class(kin.obj) == "matrix"){
      new.pheno <- kin.dat$corrected.pheno
      new.geno <- kin.dat$corrected.geno
      new.covar <- kin.dat$corrected.covar
      err.cov <- kin.dat$err.cov
    #}else{
      #taking out the LTCO option for now, since it has weird behavior
     # new.pheno <- kin.dat[[1]]$corrected.pheno
     # new.geno <- kin.dat[[1]]$corrected.geno
     # new.covar <- kin.dat[[1]]$corrected.covar
     # err.cov <- kin.dat[[1]]$err.cov
    #  marker.chr <- get.marker.chr(data.obj, m)
    #  non.covar <- setdiff(marker.chr, 0)
    #  if(length(non.covar) == 0){kin.name = "overall"}#if both markers are covariates
    #  if(length(non.covar) == 1){kin.name = paste(rep(non.covar, 2), collapse = ",")}
    #  if(length(non.covar) == 2){kin.name = paste(marker.chr, collapse = ",")}
    #  kin.locale <- which(names(kin.obj) == kin.name)
      
    #  new.pheno <- kin.dat[[kin.locale]]$corrected.pheno
    #  new.geno <- kin.dat[[kin.locale]]$corrected.geno
    #  new.covar <- kin.dat[[kin.locale]]$corrected.covar
    #  err.cov <- kin.dat[[kin.locale]]$err.cov			
    #}
    
    int.term = matrix(solve(err.cov) %*% new.geno[,m[1]]*new.geno[,m[2]], ncol = 1)
    pairscan.results <- one.pairscan.parallel(data.obj, phenotype.vector = new.pheno,
    		genotype.matrix = new.geno, int = int.term, covar.vector = new.covar, 
    		paired.markers = matrix(m, ncol = 2), n.perm = 0, verbose = verbose, 
    		run.parallel = run.parallel, n.cores = n.cores)

    if(is.null(pairscan.results[[1]])){
#      marker.num <- get.marker.num(data.obj, m)
		marker.num <- m
      dummyV <- c(marker.num, rep(NA, 3))
      results <- list("effects" = dummyV, "se" = dummyV, "cov" = c(dummyV, rep(NA, 4)))
    }else{
      results <- list("effects" = pairscan.results[[1]]$pairscan.effects, "se" = pairscan.results[[1]]$pairscan.se, 
                      "cov" = pairscan.results[[1]]$model.covariance)
    }
    return(results)
  }
  
  
  get.covar.stats <- function(m, kin.dat){
    
    #if(class(kin.obj) == "matrix"){
      new.pheno <- kin.dat$corrected.pheno
      new.geno <- kin.dat$corrected.geno
      new.covar <- kin.dat$corrected.covar
      err.cov <- kin.dat$err.cov
      #taking out LTCO option for now, since it has werid results
    #}else{
     # new.pheno <- kin.dat[[1]]$corrected.pheno
     # new.geno <- kin.dat[[1]]$corrected.geno
     # new.covar <- kin.dat[[1]]$corrected.covar
     # err.cov <- kin.dat[[1]]$err.cov
    #  marker.chr <- get.marker.chr(data.obj, m)
    #  non.covar <- setdiff(marker.chr, 0)
    #  if(length(non.covar) == 0){kin.name = "overall"}#if both markers are covariates
    #  if(length(non.covar) == 1){kin.name = paste(rep(non.covar, 2), collapse = ",")}
    #  if(length(non.covar) == 2){kin.name = paste(marker.chr, collapse = ",")}
    #  kin.locale <- which(names(kin.obj) == kin.name)
      
    #  new.pheno <- kin.dat[[kin.locale]]$corrected.pheno
    #  new.geno <- kin.dat[[kin.locale]]$corrected.geno
    #  new.covar <- kin.dat[[kin.locale]]$corrected.covar
    #  err.cov <- kin.dat[[kin.locale]]$err.cov			
    #}
    
    
    covar.locale <- which(covar.names %in% m)
    if(length(covar.locale) > 0){
      new.covar <- new.covar[,-covar.locale,drop=FALSE]
    }
    int.term = solve(err.cov) %*% new.geno[,m[1]]*new.geno[,m[2]]
    
    pairscan.results <- one.pairscan.parallel(data.obj, phenotype.vector = new.pheno, 
    		genotype.matrix = new.geno, int = int.term, 
    		covar.vector = new.covar, paired.markers = matrix(m, ncol = 2), 
        n.perm = 0, verbose = verbose, run.parallel = run.parallel, 
        n.cores = n.cores)
    
    if(is.null(pairscan.results[[1]])){
      marker.num <- get.marker.num(data.obj, m)
      dummyV <- c(marker.num, rep(NA, 3))
      results <- list("effects" = dummyV, "se" = dummyV, "cov" = c(dummyV, rep(NA, 4)))
    }else{
      results <- list("effects" = pairscan.results[[1]]$pairscan.effects, "se" = pairscan.results[[1]]$pairscan.se, "cov" = pairscan.results[[1]]$model.covariance)
    }
    return(results)
  }
  #============================================================================
  
  
  for(p in 1:num.pheno){
    if(verbose){
      cat("\nScanning phenotype ", colnames(pheno)[p], ":\n", sep = "")
    }
    
    covar.vector <- covar.info$covar.table
    pheno.vector <- pheno[,p,drop=FALSE]
    
    #sink the warnings from regress about solutions close to zero to a file
    sink(file.path(data.obj$results_path,"regress.warnings"))
    #if(class(kin.obj) == "matrix"){
      #if we are using an overall kinship matrix
      #convert to epistatic kinship matrix
      if(class(kin.obj) == "matrix"){
        kin.mat <- kin.obj
      }else{
        epi.locale <- which(names(kin.obj) == "emat")
        if(length(epi.locale) > 0){
          kin.mat <- kin.obj[[epi.locale]]
        }else{
          kin.mat <- kin.obj[[1]]
        }
      }
      
      kin.dat <- kinship.on.the.fly(kin.mat, geno, chr1 = NULL, chr2 = NULL, 
      phenoV = pheno.vector, covarV = covar.vector)
    #}else{
      #If we are using LTCO (taking this out for now)
    #  chr.pairs <- Reduce("rbind", strsplit(names(kin.obj), ","))
    #  kin.dat <- apply(chr.pairs, 1, function(x) kinship.on.the.fly(kin.obj, geno, 
    #  x[1], x[2], phenoV = pheno.vector, covarV = covar.vector))
    #}
    sink(NULL) #stop sinking output
    
    if (run.parallel) {
      cl <- parallel::makeCluster(n.cores)
      doParallel::registerDoParallel(cl)
      cape.dir.full <- find.package("cape")
      cape.dir <- str_replace(cape.dir.full,"cape_pkg/cape","cape_pkg")
      parallel::clusterExport(cl, varlist="cape.dir", envir=environment())
      parallel::clusterEvalQ(cl, .libPaths(cape.dir))
      pairscan.results <- foreach::foreach(m = t(marker.pairs), .export=c("rankMatrix", "one.pairscan.parallel", "get.covar", "get.marker.num","get.marker.chr"), .packages = 'cape') %dopar% {
                                        get.marker.pair.stats(m, kin.dat)
                                      }
      parallel::stopCluster(cl)
      
    } else {
      
     pairscan.results <- lapply(1:nrow(marker.pairs), 
     function(x) get.marker.pair.stats(m = marker.pairs[x,], kin.dat))      
      
    }
    
    effects.mat <- matrix(unlist(lapply(pairscan.results, function(x) x$effects)), nrow = nrow(marker.pairs), byrow = TRUE)
    colnames(effects.mat) <- c("marker1", "marker2", "marker1", "marker2", "marker1:marker2")
    
    se.mat <- matrix(unlist(lapply(pairscan.results, function(x) x$se)), nrow = nrow(marker.pairs), byrow = TRUE)
    colnames(se.mat) <- c("marker1", "marker2", "marker1", "marker2", "marker1:marker2")
    
    cov.mat <- matrix(unlist(lapply(pairscan.results, function(x) x$cov)), nrow = nrow(marker.pairs), byrow = TRUE)
    
    #if there are covariates to test explicitly	
    if(num.p.covar > 0){
      for(cv in 1:num.p.covar){
        if(verbose){cat("\tCovariate:", covar.names[p.covar.locale[cv]], "\n")}
        cv.marker.locale <- c(which(marker.pairs[,1] == covar.names[p.covar.locale[cv]]), which(marker.pairs[,2] == covar.names[p.covar.locale[cv]]))
        cv.markers <- marker.pairs[cv.marker.locale,,drop=FALSE]
        num.cv.pairs <- dim(cv.markers)[1]
        
        if(num.cv.pairs > 0){
          
          if (run.parallel) {
            
            cl <- parallel::makeCluster(n.cores)
            doParallel::registerDoParallel(cl)
            cape.dir.full <- find.package("cape")
            cape.dir <- str_replace(cape.dir.full,"cape_pkg/cape","cape_pkg")
            parallel::clusterExport(cl, varlist=c("rankMatrix", "one.pairscan.parallel", "get.covar", "get.marker.num", "get.marker.chr","cape.dir"), envir=environment())
            parallel::clusterEvalQ(cl, .libPaths(cape.dir))
            #.export = c("rankMatrix", "one.pairscan.parallel", "get.covar", "get.marker.num", "get.marker.chr")
            covar.results <- foreach::foreach(m = t(cv.markers), .packages = 'cape'
            		) %dopar% {
              get.covar.stats(m, kin.dat)
            }
            parallel::stopCluster(cl)
            
          } else {
            
            covar.results <- vector(mode = "list", length = num.cv.pairs)
            for (ind in 1:num.cv.pairs) {
              covar.results[[ind]] <- get.covar.stats(m = cv.markers[ind,], kin.dat)
            }
            
          }
          
          covar.effects <- matrix(unlist(lapply(covar.results, function(x) x$effects)), nrow = num.cv.pairs, byrow = TRUE)
          colnames(covar.effects) <- c("marker1", "marker2", "marker1", "marker2", "marker1:marker2")
          
          covar.se <- matrix(unlist(lapply(covar.results, function(x) x$se)), nrow = num.cv.pairs, byrow = TRUE)
          colnames(covar.se) <- c("marker1", "marker2", "marker1", "marker2", "marker1:marker2")
          covar.cov <- matrix(unlist(lapply(covar.results, function(x) x$cov)), nrow = num.cv.pairs, byrow = TRUE)
          
          effects.mat <- rbind(effects.mat, covar.effects)
          se.mat <- rbind(se.mat, covar.se)
          cov.mat <- rbind(cov.mat, covar.cov)
          
        }#end case for when there are pairs of markers with covariates
      } #end looping through markers paired with covariates
    }#end case for when there are phenotypic covariates to test
    
  
    
    
    #add the results for the phenotype
    pheno.results <- list(effects.mat, se.mat, cov.mat)
    names(pheno.results) <- c("pairscan.effects", "pairscan.se", "model.covariance")
    results.list[[p]] <- pheno.results
  }	#end looping over phenotypes	
  
  unlink(file.path(data.obj$results_path,"regress.warnings"))
  return(results.list)
  
  
}