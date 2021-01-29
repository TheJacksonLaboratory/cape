#' Run the pairscan with a kinship correction
#' 
#' This function is called by \code{\link{pairscan}}
#' when a kinship correction is requested. It adjusts 
#' each variable according to the kinship matrix using
#' \code{\link{kin_adjust}} and then fits linear
#' pairwise models to the adjusted data.
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' @param scan_what A character string uniquely identifying whether eigentraits
#'   or raw traits should be scanned. Options are "eigentraits", "raw_traits"
#' @param marker_pairs A two-column matrix containing the marker pairs
#' to be tested.
#' @param kin_obj a kinship object
#' @param verbose A logical value indicating whether to 
#' print progress to the screen
#' @param run_parallel A logical value indicating 
#' whether parallel processing should be used
#' @param n_cores The number of cores to be used if run_parallel is TRUE
#' 
#' @return This function returns a list with three elements. 
#' The elements contain the marker pair effect sizes, the marker
#' pair standard errors, and the covariance matrix for each test.
#' The output is then further processed by \code{\link{pairscan}}.
#' 
#' @import parallel
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @keywords internal
#' 
#' 
pairscan_kin <- function(data_obj, geno_obj, scan_what, marker_pairs, 
kin_obj, verbose = FALSE, run_parallel = FALSE, n_cores = 2){
  
  m = NULL #for appeasing R CMD check
  
  pheno <- get_pheno(data_obj, scan_what)	
  geno <- get_geno_with_covar(data_obj, geno_obj, for_pairscan = TRUE)
  num_pheno <- dim(pheno)[2]
  results_list <- vector(mode = "list", length = num_pheno)
  names(results_list) <- colnames(pheno)
  
  covar_info <- get_covar(data_obj)
  is_char <- as.logical(is.na(suppressWarnings(as.numeric(marker_pairs[1,1]))))
  if(is_char){
    covar_names <- get_marker_name(data_obj, covar_info$covar_names)
  }else{
    covar_names <- get_marker_num(data_obj, covar_info$covar_names)	
  }
  
  p_covar_locale <- which(covar_info$covar_type == "p")
  num_p_covar <- length(p_covar_locale)
  #============================================================================
  #internal functions
  #============================================================================
  
  get_marker_pair_stats <- function(m, kin_dat){
    
    if(data_obj$kinship_type == "overall"){
      new_pheno <- kin_dat$corrected_pheno
      new_geno <- kin_dat$corrected_geno
      new_covar <- kin_dat$corrected_covar
      err_cov <- kin_dat$err_cov
    }else{
      marker_chr <- get_marker_chr(data_obj, m)      
      non_covar <- setdiff(marker_chr, 0)

      if(length(non_covar) == 0){kin_name = "overall"}#if both markers are covariates
      if(length(non_covar) == 1){kin_name = paste(rep(non_covar, 2), collapse = ",")}
      if(length(non_covar) == 2){kin_name = paste(marker_chr, collapse = ",")}
      
      kin_locale <- which(names(kin_obj) == kin_name)
      
      new_pheno <- kin_dat[[kin_locale]]$corrected_pheno
      new_geno <- kin_dat[[kin_locale]]$corrected_geno
      new_covar <- kin_dat[[kin_locale]]$corrected_covar
      err_cov <- kin_dat[[kin_locale]]$err_cov			
    }
    
    int_term = matrix(solve(err_cov) %*% new_geno[,m[1]]*new_geno[,m[2]], ncol = 1)
    pairscan_results <- one_pairscan_parallel(data_obj, phenotype_vector = new_pheno,
    		genotype_matrix = new_geno, int = int_term, covar_vector = new_covar, 
    		paired_markers = matrix(m, ncol = 2), n_perm = 0, verbose = verbose, 
    		run_parallel = run_parallel, n_cores = n_cores)

    if(is.null(pairscan_results[[1]])){
#      marker_num <- get_marker_num(data_obj, m)
		marker_num <- m
      dummyV <- c(marker_num, rep(NA, 3))
      results <- list("effects" = dummyV, "se" = dummyV, "cov" = c(dummyV, rep(NA, 4)))
    }else{
      results <- list("effects" = pairscan_results[[1]]$pairscan_effects, "se" = pairscan_results[[1]]$pairscan_se, 
                      "cov" = pairscan_results[[1]]$model_covariance)
    }
    return(results)
  }
  
  
  get_covar_stats <- function(m, kin_dat){
    
    if(data_obj$kinship_type == "overall"){
      new_pheno <- kin_dat$corrected_pheno
      new_geno <- kin_dat$corrected_geno
      new_covar <- kin_dat$corrected_covar
      err_cov <- kin_dat$err_cov
    }else{
      marker_chr <- get_marker_chr(data_obj, m)
      non_covar <- setdiff(marker_chr, 0)
      if(length(non_covar) == 0){kin_name = "overall"}#if both markers are covariates
      if(length(non_covar) == 1){kin_name = paste(rep(non_covar, 2), collapse = ",")}
      if(length(non_covar) == 2){kin_name = paste(marker_chr, collapse = ",")}
      kin_locale <- which(names(kin_obj) == kin_name)
      
      new_pheno <- kin_dat[[kin_locale]]$corrected_pheno
      new_geno <- kin_dat[[kin_locale]]$corrected_geno
      new_covar <- kin_dat[[kin_locale]]$corrected_covar
      err_cov <- kin_dat[[kin_locale]]$err_cov			
    }
    
    covar_locale <- which(covar_names %in% m)
    if(length(covar_locale) > 0){
      new_covar <- new_covar[,-covar_locale,drop=FALSE]
    }
    int_term = solve(err_cov) %*% new_geno[,m[1]]*new_geno[,m[2]]
    
    pairscan_results <- one_pairscan_parallel(data_obj, phenotype_vector = new_pheno, 
    		genotype_matrix = new_geno, int = int_term, 
    		covar_vector = new_covar, paired_markers = matrix(m, ncol = 2), 
        n_perm = 0, verbose = verbose, run_parallel = run_parallel, 
        n_cores = n_cores)
    
    if(is.null(pairscan_results[[1]])){
      marker_num <- get_marker_num(data_obj, m)
      dummyV <- c(marker_num, rep(NA, 3))
      results <- list("effects" = dummyV, "se" = dummyV, "cov" = c(dummyV, rep(NA, 4)))
    }else{
      results <- list("effects" = pairscan_results[[1]]$pairscan_effects, "se" = pairscan_results[[1]]$pairscan_se, "cov" = pairscan_results[[1]]$model_covariance)
    }
    return(results)
  }
  #============================================================================
  
  
  for(p in 1:num_pheno){
    if(verbose){
      cat("\nScanning phenotype ", colnames(pheno)[p], ":\n", sep = "")
    }
    
    covar_vector <- covar_info$covar_table
    pheno_vector <- pheno[,p,drop=FALSE]
    
    #sink the warnings from regress about solutions close to zero to a file
    sink(file.path(data_obj$results_path,"regress.warnings.pairscan"))
    
    if(data_obj$kinship_type == "overall"){
      #if we are using an overall kinship matrix
        kin_mat <- kin_obj
        kin_dat <- kin_adjust(kin_mat, geno, chr1 = NULL, chr2 = NULL, 
        phenoV = pheno_vector, covarV = covar_vector)
    }else{
      #If we are using LTCO
      chr_pairs <- Reduce("rbind", strsplit(names(kin_obj), ","))
      kin_dat <- apply(chr_pairs, 1, function(x) kin_adjust(kin_obj, geno, 
      x[1], x[2], phenoV = pheno_vector, covarV = covar_vector))
      names(kin_dat) <- names(kin_obj)
    }
    sink(NULL) #stop sinking output
    
    if (run_parallel) {
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      cape_dir_full <- find.package("cape")
      cape_dir <- str_replace(cape_dir_full,"cape_pkg/cape","cape_pkg")
      clusterExport(cl, varlist="cape_dir", envir=environment())
      clusterEvalQ(cl, .libPaths(cape_dir))
      pairscan_results <- foreach(m = t(marker_pairs), .export=c("rankMatrix", "one_pairscan_parallel", "get_covar", "get_marker_num","get_marker_chr"), .packages = 'cape') %dopar% {
                                        get_marker_pair_stats(m, kin_dat)
                                      }
      stopCluster(cl)
      
    } else {
      
     pairscan_results <- lapply(1:nrow(marker_pairs), 
     function(x) get_marker_pair_stats(m = marker_pairs[x,], kin_dat))      
      
    }
    
    effects_mat <- matrix(unlist(lapply(pairscan_results, function(x) x$effects)), nrow = nrow(marker_pairs), byrow = TRUE)
    colnames(effects_mat) <- c("marker1", "marker2", "marker1", "marker2", "marker1:marker2")
    
    se_mat <- matrix(unlist(lapply(pairscan_results, function(x) x$se)), nrow = nrow(marker_pairs), byrow = TRUE)
    colnames(se_mat) <- c("marker1", "marker2", "marker1", "marker2", "marker1:marker2")
    
    cov_mat <- matrix(unlist(lapply(pairscan_results, function(x) x$cov)), nrow = nrow(marker_pairs), byrow = TRUE)
    
    #if there are covariates to test explicitly	
    if(num_p_covar > 0){
      for(cv in 1:num_p_covar){
        if(verbose){cat("\tCovariate:", covar_names[p_covar_locale[cv]], "\n")}
        cv_marker_locale <- c(which(marker_pairs[,1] == covar_names[p_covar_locale[cv]]), which(marker_pairs[,2] == covar_names[p_covar_locale[cv]]))
        cv_markers <- marker_pairs[cv_marker_locale,,drop=FALSE]
        num_cv_pairs <- dim(cv_markers)[1]
        
        if(num_cv_pairs > 0){
          
          if (run_parallel) {
            
            cl <- makeCluster(n_cores)
            registerDoParallel(cl)
            cape_dir_full <- find.package("cape")
            cape_dir <- str_replace(cape_dir_full,"cape_pkg/cape","cape_pkg")
            clusterExport(cl, varlist=c("rankMatrix", "one_pairscan_parallel", "get_covar", "get_marker_num", "get_marker_chr","cape_dir"), envir=environment())
            clusterEvalQ(cl, .libPaths(cape_dir))
            #.export = c("rankMatrix", "one_pairscan_parallel", "get_covar", "get_marker_num", "get_marker_chr")
            covar_results <- foreach(m = t(cv_markers), .packages = 'cape'
            		) %dopar% {
              get_covar_stats(m, kin_dat)
            }
            stopCluster(cl)
            
          } else {
            
            covar_results <- vector(mode = "list", length = num_cv_pairs)
            for (ind in 1:num_cv_pairs) {
              covar_results[[ind]] <- get_covar_stats(m = cv_markers[ind,], kin_dat)
            }
            
          }
          
          covar_effects <- matrix(unlist(lapply(covar_results, function(x) x$effects)), nrow = num_cv_pairs, byrow = TRUE)
          colnames(covar_effects) <- c("marker1", "marker2", "marker1", "marker2", "marker1:marker2")
          
          covar_se <- matrix(unlist(lapply(covar_results, function(x) x$se)), nrow = num_cv_pairs, byrow = TRUE)
          colnames(covar_se) <- c("marker1", "marker2", "marker1", "marker2", "marker1:marker2")
          covar_cov <- matrix(unlist(lapply(covar_results, function(x) x$cov)), nrow = num_cv_pairs, byrow = TRUE)
          
          effects_mat <- rbind(effects_mat, covar_effects)
          se_mat <- rbind(se_mat, covar_se)
          cov_mat <- rbind(cov_mat, covar_cov)
          
        }#end case for when there are pairs of markers with covariates
      } #end looping through markers paired with covariates
    }#end case for when there are phenotypic covariates to test
    
  
    
    
    #add the results for the phenotype
    pheno_results <- list(effects_mat, se_mat, cov_mat)
    names(pheno_results) <- c("pairscan_effects", "pairscan_se", "model_covariance")
    results_list[[p]] <- pheno_results
  }	#end looping over phenotypes	
  
  #unlink(file_path(data_obj$results_path,"regress_warnings"))
  return(results_list)
  
  
}