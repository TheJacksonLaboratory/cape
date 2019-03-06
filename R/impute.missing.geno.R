#' impute missing genotype data using k nearest neighbors
#' 
#' This function uses k nearest neighbors to impute missing genotype data
#' on a per chromosome basis prioritize lets you remove individuals before
#' markers, vice versa, or decide based on which will remove fewer elements.
#' If you have a separate data.obj and geno.obj, this script returns both 
#' as \code{list(data.obj, geno.obj)}. The data.obj has the updated markers and 
#' individuals, and the geno.obj has the imputed genetic data. These objects
#' must be separated after running this script. max.region.size specifies the
#' maximum number of markers to be used in calculating individual similarity.
#' There is a trade-off between the time it takes to calculate a distance
#' matrix for a large matrix and the time it takes to slide through the
#' genome imputing markers. This function does not yet support imputation of
#' covariates.
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param geno.obj a genotype object
#' @param k integer number of neighboring markers, default = 10
#' @param ind.missing.thresh percent of individuals that are acceptible to remove, default = 0
#' @param marker.missing.thresh percent of genotype markers that are acceptible to remove, default = 0
#' @param prioritize the basis prioritization is one of 
#'        "fewer" = remove the fewest possible cells from the matrix
#'        "ind" = remove the fewest possible individuals
#'        "marker" = remove the fewest possible markers
#' @param max.region.size maximum number of markers to be used in calculating individual similarity
#' @param min.region.size minimum number of markers to be used in calculating individual similarity
#' @param run.parallel boolean
#' @param n.cores integer number of available CPU cores to use for parallel processing
#' @param verbose boolean
#'
#' @return an updated cape object
#'
#' @export
impute.missing.geno <- function(data.obj, geno.obj = NULL, k = 10, ind.missing.thresh = 0, 
                                marker.missing.thresh = 0, prioritize = c("fewer", "ind", "marker"),
                                max.region.size = NULL, min.region.size = NULL, run.parallel = TRUE,
                                verbose = FALSE, n.cores = 2){
  
  impute.full.genome = FALSE
  p.choice <- grep("f", prioritize)
  
  if(length(p.choice) > 0){prioritize <- "fewer"}
  
  if(!impute.full.genome){
    new.geno <- get.geno(data.obj, geno.obj)
    chr <- data.obj$chromosome
    u_chr <- unique(chr)
    ind.used <- rownames(data.obj$pheno)
    ind.used.locale <- match(ind.used, rownames(new.geno))
    new.geno <- new.geno[ind.used.locale,,]
  }else{
    new.geno <- geno.obj$geno	
    if(is.null(new.geno)){
      stop("I could not find a genotype matrix in the genotype object. Please check or set impute.full.genome to FALSE")
    }
    chr <- data.obj$chromosome
    u_chr <- unique(chr)
    ind.used <- rownames(data.obj$pheno)
    ind.used.locale <- match(ind.used, rownames(new.geno))
    new.geno <- new.geno[ind.used.locale,,]
  }
  
  chr.lengths <- apply(matrix(u_chr, nrow = 1), 2, function(x) length(which(chr == x)))
  
  #go through the chromosome at the window size specified and impute genotypes
  if(is.null(max.region.size) || max.region.size > max(chr.lengths)){
    max.region <- max(chr.lengths)
    min.region <- min(chr.lengths)
  }else{
    max.region <- max.region.size
    min.region <- min.region.size
    if(is.null(min.region)){
      min.region <- 10			
    }
  }
  if(max.region > 1000){
    choice <- 0
    while(choice != "y" && choice != "n"){
      choice <- readline(prompt = "At least one chromosome has more than 1000 markers.\n\tCalculation of large distance matrices may crash R.\n\tSet max.region.size to specify maximum number of markers in imputed region.\n\tWould you like to continue with the current settings? (y/n)")
      if(choice == "n"){stop()}
    }	
  }
  
  
  marker.names <- dimnames(new.geno)[[3]]
  ind.names <- rownames(new.geno)
  
  num.sections <- length(seq(1, dim(new.geno)[2], max.region))
  geno.chunks <- vector(mode = "list", length = num.sections)
  
  if(verbose){cat("Dividing genome into chunks.\n")}
  chunk.pos <- 1
  for(ch in u_chr){
    if(verbose){report.progress(which(u_chr == ch), length(u_chr))}
    #pull out all the genotypes for the chromosome
    chr.locale <- which(chr == ch)
    if(length(chr.locale) <= max.region){
      geno.chunks[[chunk.pos]] <- new.geno[,,chr.locale,drop=FALSE]
      chunk.pos <- chunk.pos + 1
    }else{
      marker.seq <- seq(1, length(chr.locale), max.region)
      if(tail(marker.seq, 1) < length(chr.locale)){
        if(length(chr.locale)- tail(marker.seq,1) >= min.region){
          marker.seq <- c(marker.seq, length(chr.locale))
        }else{
          marker.seq[length(marker.seq)] <- length(chr.locale)
        }
      }
      for(m in 1:(length(marker.seq)-1)){
        if(m == length(marker.seq)-1){
          last.pos <- chr.locale[(marker.seq[m+1])]
        }else{
          last.pos <- chr.locale[(marker.seq[m+1]-1)]	
        }
        geno.chunks[[chunk.pos]] <- new.geno[,,chr.locale[marker.seq[m]]:last.pos]
        chunk.pos <- chunk.pos + 1
      }
    }
  }
  if(verbose){cat("\n")}
  
  impute.section <- function(sec.geno){
    if(dim(sec.geno)[2] == 1){return(sec.geno)}
    ind.missing.geno <- which(apply(sec.geno, 1, function(x) length(which(is.na(x)))) > 0)
    ind.dist <- array(NA, dim = c(nrow(sec.geno), nrow(sec.geno), dim(sec.geno)[2]))
    if(length(ind.missing.geno) > 0){
      for(l in 1:dim(sec.geno)[2]){
        ind.dist[,,l] <- as.matrix(dist(sec.geno[,l,]))
      }
      ind.dist.mean <- flatten.array(ind.dist, 1, 2, function(x) mean(x, na.rm = TRUE))
      dimnames(ind.dist.mean) <- list(dimnames(sec.geno)[[1]], dimnames(sec.geno)[[1]])
      for(id in 1:length(ind.missing.geno)){
        ind <- sec.geno[ind.missing.geno[id],,,drop=FALSE]
        neighbors <- sort(ind.dist.mean[ind.missing.geno[id],])
        if(length(neighbors) < 2){
          next()
        }
        nearest.neighbors <- names(neighbors[2:min(c((k+1), length(neighbors)))])
        neighbor.weights <- neighbors[nearest.neighbors]/sum(neighbors[nearest.neighbors], na.rm = TRUE)
        missing.geno <- which(is.na(ind), arr.ind = TRUE)
        
        
        missing.geno.ind <- array(sec.geno[as.character(nearest.neighbors),missing.geno[,2],missing.geno[,3],drop=FALSE], dim = c(length(nearest.neighbors), ncol(sec.geno), length(unique(missing.geno[,3]))))
        
        
        weight.mat <- array(neighbor.weights, dim = dim(missing.geno.ind))
        missing.geno.weighted <- missing.geno.ind*weight.mat
        ind[missing.geno] <- matrix(apply(missing.geno.weighted, 3, function(x) colSums(x, na.rm = TRUE)), nrow = dim(ind)[2], ncol = dim(missing.geno.ind)[3])
        sec.geno[ind.missing.geno[id],,] <- ind
      } #end looping through individuals with missing genotypes in this section
    } #end case for if there are individuals with missing genotypes
    return(sec.geno)
  } #end impute.section
  
  browser()
  
  if(run.parallel){
    if(verbose){cat("Imputing missing genotypes...\n")}
    cl <- makeCluster(n.cores)
    registerDoParallel(cl)
    imputed.geno <- foreach(m = geno.chunks, .export = "flatten.array") %dopar% {
      impute.section(m)
    }
    stopCluster(cl)
  }else{
    # if(verbose){cat("Imputing missing genotypes...\n")
    # imputed.geno <- lapply_pb(geno.chunks, impute.section)
    # }else{-
    imputed.geno <- lapply(geno.chunks, impute.section)
    # }
  }
  
  if(verbose){cat("Rebuilding genotype matrix...\n")}
  #now unlist the result into the new genotype matrix
  new.col <- sum(unlist(lapply(imputed.geno, function(x) dim(x)[3])))
  imp.geno <- array(unlist(imputed.geno), dim = c(nrow(new.geno), ncol(new.geno), new.col))
  dimnames(imp.geno) <- dimnames(new.geno)
  
  if(verbose){cat("Removing missing data...\n")}
  new.geno <- remove.missing.genotype.data(data.obj, imp.geno, ind.missing.thresh = 0, marker.missing.thresh = 0, prioritize = c("ind", "marker", "fewer"))
  
  
  final.obj <- list("data.obj" = new.geno, "geno.obj" = imp.geno)
  return(final.obj)
  
}