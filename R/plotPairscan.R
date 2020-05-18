#' Plot the result of the 2D scan
#'
#' This script plots the results of the 2D scan.
#' It plots a 2D matrix of the t-statistics of
#' the interactions between all pairs of 
#' markers.
#' The data object is only used to get marker
#' order. This might be redundant, but it seemed
#' easier for the user to put this in then to 
#' extract the marker order themselves
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param pairscan.obj a pariscan object from \code{\link{pairscan}}
#' @param phenotype default = NULL, can be drawn from \code{pairscan.obj$pairscan.results}
#' @param standardized default = FALSE
#' @param show.marker.labels
#' @param show.chr
#' @param label.chr
#' @param show.alleles
#' @param allele.labels
#' @param pos.col default = "brown"
#' @param neg.col default = "blue"
#' @param color.scheme options are "DO/CC" and "other"
#' @param pdf.label default = "Pairscan.Regression.pdf"
#'
plotPairscan <- function(data.obj, pairscan.obj, phenotype = NULL, standardized = FALSE, show.marker.labels = FALSE,
                         show.chr = TRUE, label.chr = TRUE, show.alleles = TRUE, allele.labels = NULL, pos.col = "brown",
                         neg.col = "blue", color.scheme = c("DO/CC", "other"), pdf.label = "Pairscan.Regression.pdf") {
  
  require(igraph)
  
  pairscan.results <- pairscan.obj$pairscan.results
  
  marker.pairs <- pairscan.results[[1]][[1]][,1:2]
  #get the markers used in the pair scan and sort them.
  markers <- unique(as.vector(marker.pairs))
  markers.no.allele <- unlist(lapply(strsplit(markers, "_"), function(x) x[[1]][1]))
  
  sorted.markers <- markers[order(as.numeric(get.marker.num(data.obj, markers.no.allele)))]
  
  #get coordinates of the chromosome boundaries
  if(show.chr){
    chromosomes <- get.marker.chr(data.obj, markers.no.allele)
    u_chr <- unique(chromosomes[!is.na(chromosomes)])
    chr.boundaries <- apply(matrix(u_chr, ncol = 1), 1, function(x) max(which(chromosomes == x))) + 0.5
    chr.boundaries <- c(0, chr.boundaries)
    if(label.chr){
      chr.names <- unique(chromosomes)
    }else{
      chr.names <- NULL	
    }
  }else{
    chr.boundaries <- NULL
    chr.names <- NULL
  }
  
  
  if(show.alleles){
    
    alleles <- unique(sapply(strsplit(colnames(data.obj$geno_for_pairscan), "_"), function(x) x[2]))
    allele.colors <- get.allele.colors(color.scheme, alleles)
    allele.labels <- allele.colors[,2]
    
    
    all.alleles <- unlist(lapply(strsplit(sorted.markers, "_"), function(x) x[2]))
    allele.cols <- allele.colors[match(all.alleles, alleles),3]
    
  }else{
    allele.cols <- NULL
  }
  
  pairscan.result <- pairscan.obj$pairscan.results
  
  if(is.null(pairscan.result)){
    stop("pairscan() must be run before plotPairscan()")
  }
  
  if(is.null(phenotype)){
    phenotype <- names(pairscan.result)
  }
  
  pheno.num <- which(names(pairscan.result) %in% phenotype)
  
  if(length(pheno.num) < length(phenotype)){
    not.found <- which(!(phenotype %in% names(pairscan.result)))
    cat("I couldn't find the following phenotypes:")
    cat(phenotype[not.found], sep = "\n")
    stop()
  }
  
  #collect the results, so we can put them on the same scale
  all.results.mats <- list()
  min.x <- 0
  max.x <- 0
  #for each phenotype scanned
  for(p in pheno.num){
    #build a results matrix
    results.mat <- matrix(0, length(markers), length(markers))
    colnames(results.mat) <- rownames(results.mat) <- sorted.markers
    
    if(standardized){
      pair.int <- as.numeric(pairscan.result[[p]][[1]][, "marker1:marker2"])/as.numeric(pairscan.result[[p]][[2]][,"marker1:marker2"])
    }else{
      pair.int <- as.numeric(pairscan.result[[p]][[1]][,"marker1:marker2"])
    }
    
    #create symetric matrices with the values
    net <- graph.edgelist(pairscan.result[[p]][[1]][,1:2])
    E(net)$weight <- pair.int
    results.mat <- as.matrix(as_adjacency_matrix(net, attr = "weight"))
    results.mat <- results.mat + t(results.mat)
    
    diag(results.mat) <- NA
    all.results.mats[[p]] <- results.mat
    min.x <- max(abs(c(max.x, results.mat)), na.rm = TRUE)*-1
    max.x <- max(abs(c(max.x, results.mat)), na.rm = TRUE)
    
  }	
  
  if (endsWith(pdf.label, "pdf")) {
    pdf(pdf.label, width = 7, height = 6)
  } else if (endsWith(pdf.label, "jpg")) {
    jpeg(pdf.label, quality = 100)
  }
  
  for(p in 1:length(pheno.num)){
    myImagePlot(x = all.results.mats[[pheno.num[p]]], xlab = "marker1", ylab = "marker2", main = phenotype[p], min.x = min.x, max.x = max.x, show.labels = show.marker.labels, chromosome.coordinates = chr.boundaries, chr.names = chr.names, allele.cols = allele.cols, pos.col = pos.col, neg.col = neg.col)
  }
  dev.off()
  
}