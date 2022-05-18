#' Impute missing genotype data using k nearest neighbors
#' 
#' This function uses k nearest neighbors to impute missing genotype data
#' on a per chromosome basis. If missing genotypes remain after imputations
#' the user can prioritize whether to remove individuals, markers, or whichever
#' has fewer missing values.
#' 
#' This function is run by \code{\link{run_cape}} and runs automatically if
#' a kinship correction is specified and there are missing values in the 
#' genotype object.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' @param k The number of nearest neighbors to use to impute missing data. Defaults to 10.
#' @param ind_missing_thresh percent A percentage of acceptable missing data. After imputation
#' if an individual is missing more data than the percent specified, it will be removed.
#' @param marker_missing_thresh A percentage of acceptable missing data. After imputation
#' if a marker is missing more data than the percent specified, it will be removed.
#' @param prioritize How to prioritize removal of rows and columns with missing data.
#' "ind" = remove individuals with missing data exceeding the threshold before 
#' considering markers to remove.
#' "marker" = remove markers with missing data exceeding the threshold before
#' considering individuals to remove.
#' "fewer" = Determine how much data will be removed by prioritizing individuals
#' or markers. Remove data in whichever order removes the least amount of data.
#' @param max_region_size maximum number of markers to be used in calculating individual similarity.
#' Defaults to the minimum chromosome size.
#' @param min_region_size minimum number of markers to be used in calculating individual similarity
#' Defaults to the maximum chromosome size.
#' @param run_parallel A logical value indicating whether to run the process in parallel
#' @param n_cores integer number of available CPU cores to use for parallel processing
#' @param verbose A logical value indicating whether to print progress to the screen.
#' 
#' @details The prioritize parameter can be a bit confusing. If after imputation,
#' there is one marker for which all data are missing, it makes sense to remove that
#' one marker rather than all individuals with missing data, since all individuals
#' would be removed. Similarly, if there is one individual with massive amounts of 
#' missing data, it makes sense to remove that individual, rather than all markers
#' that individual is missing. We recommend always using the default "fewer" option 
#' here unless you know for certain that you want to prioritize individuals or markers 
#' for removal.
#' There is no need to specify max_region_size or min_region_size, but advanced
#' users may want to specify them. There is a trade-off between the time it takes 
#' to calculate a distance matrix for a large matrix and the time it takes to slide 
#' through the genome imputing markers. This function does not yet support imputation 
#' of covariates.
#' If individuals are genotyped very densely, the user may want to specify max_region_size
#' to be smaller than the maximum chromosome size to speed calculation of similarity matrices.
#'
#' @return This function returns a list that includes both the data_obj and geno_obj
#' These objects must then be separated again to continue through the cape analysis.
#'
#' @import parallel
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' 
#' @examples 
#' \dontrun{
#' combined_obj <- impute_missing_geno(data_obj, geno_obj)
#' new_data_obj <- combined_obj$data_obj
#' noew_geno_obj <- combined_obj$geno_obj
#' }
#'
#' @export

impute_missing_geno <- function(data_obj, geno_obj = NULL, k = 10, ind_missing_thresh = 0, 
                                marker_missing_thresh = 0, prioritize = c("fewer", "ind", "marker"),
                                max_region_size = NULL, min_region_size = NULL, run_parallel = FALSE,
                                verbose = FALSE, n_cores = 2){
  
  impute_full_genome = FALSE
  p.choice <- grep("f", prioritize)
  
  if(length(p.choice) > 0){prioritize <- "fewer"}
  
  if(!impute_full_genome){
    new_geno <- get_geno(data_obj, geno_obj)
    chr <- data_obj$chromosome
    u_chr <- unique(chr)
    ind_used <- rownames(data_obj$pheno)
    ind_used_locale <- match(ind_used, rownames(new_geno))
    new_geno <- new_geno[ind_used_locale,,]
  }else{
    new_geno <- geno_obj$geno	
    if(is.null(new_geno)){
      stop("I could not find a genotype matrix in the genotype object. Please check or set impute_full_genome to FALSE")
    }
    chr <- data_obj$chromosome
    u_chr <- unique(chr)
    ind_used <- rownames(data_obj$pheno)
    ind_used_locale <- match(ind_used, rownames(new_geno))
    new_geno <- new_geno[ind_used_locale,,]
  }
  
  chr_lengths <- apply(matrix(u_chr, nrow = 1), 2, function(x) length(which(chr == x)))
  
  #go through the chromosome at the window size specified and impute genotypes
  if(is.null(max_region_size) || max_region_size > max(chr_lengths)){
    max_region <- max(chr_lengths)
    min_region <- min(chr_lengths)
  }else{
    max_region <- max_region_size
    min_region <- min_region_size
    if(is.null(min_region)){
      min_region <- 10			
    }
  }
  if(max_region > 1000){
    choice <- 0
    while(choice != "y" && choice != "n"){
      choice <- readline(prompt = "At least one chromosome has more than 1000 markers.\n\tCalculation of large distance matrices may crash R.\n\tSet max_region_size to specify maximum number of markers in imputed region.\n\tWould you like to continue with the current settings? (y/n)")
      if(choice == "n"){stop()}
    }	
  }
  
  num_sections <- length(seq(1, dim(new_geno)[2], max_region))
  geno_chunks <- vector(mode = "list", length = num_sections)
  
  if(verbose){cat("Dividing genome into chunks.\n")}
  chunk_pos <- 1
  for(ch in u_chr){
    if(verbose){report_progress(which(u_chr == ch), length(u_chr))}
    #pull out all the genotypes for the chromosome
    chr_locale <- which(chr == ch)
    if(length(chr_locale) <= max_region){
      geno_chunks[[chunk_pos]] <- new_geno[,,chr_locale,drop=FALSE]
      chunk_pos <- chunk_pos + 1
    }else{
      marker_seq <- seq(1, length(chr_locale), max_region)
      if(tail(marker_seq, 1) < length(chr_locale)){
        if(length(chr_locale)- tail(marker_seq,1) >= min_region){
          marker_seq <- c(marker_seq, length(chr_locale))
        }else{
          marker_seq[length(marker_seq)] <- length(chr_locale)
        }
      }
      for(m in 1:(length(marker_seq)-1)){
        if(m == length(marker_seq)-1){
          last_pos <- chr_locale[(marker_seq[m+1])]
        }else{
          last_pos <- chr_locale[(marker_seq[m+1]-1)]	
        }
        geno_chunks[[chunk_pos]] <- new_geno[,,chr_locale[marker_seq[m]]:last_pos]
        chunk_pos <- chunk_pos + 1
      }
    }
  }
  if(verbose){cat("\n")}
  
  impute_section <- function(sec_geno){
    if(dim(sec_geno)[2] == 1){return(sec_geno)}
    ind_missing_geno <- which(apply(sec_geno, 1, function(x) length(which(is.na(x)))) > 0)
    ind_dist <- array(NA, dim = c(nrow(sec_geno), nrow(sec_geno), dim(sec_geno)[2]))
    if(length(ind_missing_geno) > 0){
      list_geno <- lapply(1:dim(sec_geno)[2], function(x) sec_geno[,x,])
      flat_geno <- Reduce("cbind", list_geno)
      #dim(flat_geno)
      ind_dist <- as.matrix(dist(flat_geno))
      for(id in 1:length(ind_missing_geno)){
        ind <- sec_geno[ind_missing_geno[id],,,drop=FALSE]
        neighbors <- sort(ind_dist[ind_missing_geno[id],])
        if(length(neighbors) < 2){
          next()
        }
        nearest_neighbors <- names(neighbors[2:min(c((k+1), length(neighbors)))])
        neighbor_weights <- neighbors[nearest_neighbors]/sum(neighbors[nearest_neighbors], na.rm = TRUE)
        missing_geno <- which(is.na(ind), arr.ind = TRUE)
        
        missing_geno_ind <- array(sec_geno[as.character(nearest_neighbors),missing_geno[,2],missing_geno[,3],drop=FALSE], dim = c(length(nearest_neighbors), ncol(sec_geno), length(unique(missing_geno[,3]))))
        
        weight_mat <- array(neighbor_weights, dim = dim(missing_geno_ind))
        missing_geno_weighted <- missing_geno_ind*weight_mat
        ind[missing_geno] <- matrix(apply(missing_geno_weighted, 3, function(x) colSums(x, na.rm = TRUE)), nrow = dim(ind)[2], ncol = dim(missing_geno_ind)[3])
        sec_geno[ind_missing_geno[id],,] <- ind
      } #end looping through individuals with missing genotypes in this section
    } #end case for if there are individuals with missing genotypes
    return(sec_geno)
  } #end impute_section
  
  if(run_parallel){
    if(verbose){cat("Imputing missing genotypes...\n")}
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    cape_dir_full <- find.package("cape")
    cape_dir <- gsub("cape_pkg/cape", "cape_pkg", cape_dir_full)
    clusterExport(cl, "cape_dir", envir=environment())
    clusterEvalQ(cl, .libPaths(cape_dir))
    # the following line adds package variables to the parallel worker environments
    #parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
    imputed_geno <- foreach(m = geno_chunks, .export = "flatten_array") %dopar% {
      impute_section(m)
    }
    stopCluster(cl)
  }else{
    if(verbose){cat("Imputing missing genotypes...\n")}
      imputed_geno <- lapply(geno_chunks, impute_section)
  }
  
  if(verbose){cat("Rebuilding genotype matrix...\n")}
  #now unlist the result into the new genotype matrix
  new_col <- sum(unlist(lapply(imputed_geno, function(x) dim(x)[3])))
  imp_geno <- array(unlist(imputed_geno), dim = c(nrow(new_geno), ncol(new_geno), new_col))
  dimnames(imp_geno) <- dimnames(new_geno)
  
  if(verbose){cat("Removing missing data...\n")}
  new_geno <- remove_missing_genotype_data(data_obj, imp_geno, ind_missing_thresh = 0, 
  marker_missing_thresh = 0, prioritize = c("ind", "marker", "fewer"))
  
  
  final_obj <- list("data_obj" = new_geno, "geno_obj" = imp_geno)
  return(final_obj)
  
}