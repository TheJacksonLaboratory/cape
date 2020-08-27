#' Plot the result of the pairwise scan
#'
#' This function plots the results of the pairwise scan.
#' It plots a matrix of the the interactions between all 
#' pairs of markers.
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param pairscan.obj a pariscan object from \code{\link{pairscan}}
#' @param phenotype The names of the phenotypes to be plotted. If NULL,
#' all phenotypes are plotted.
#' @param standardized If TRUE, the standardized effects are plotted.
#' IF FALSE, the effect sizes are plotted.
#' @param show.marker.labels If TRUE marker labels are plotted along the
#' axes. If FALSE, they are omitted.
#' @param show.chr If TRUE, the chromosome boundaries are shown
#' @param label.chr If TRUE, the chromsomes are labeled
#' @param show.alleles If TRUE, the allele of each marker is indicated by color.
#' @param allele.labels Labels for the alleles if other than those stored in the
#' data object.
#' @param pos.col The color to use for positive main effects and interactions
#' must be one of "green", "purple", "red", "orange", "blue", "brown", "yellow", "gray"
#' see \code{\link{get.color}}
#' @param neg.col The color to use for negative main effects and interactions
#' takes the same values as pos.col.
#' @param color.scheme either "DO/CC" or "other". "DO/CC" uses the official "DO/CC"
#' colors for the DO/CC alleles  \url{https://compgen.unc.edu/wp/?page_id=577}
#' "other" uses an unrelated color pallette for multiple alleles.
#' @param pdf.label Label for the resulting file. Defaults to "Pairscan.Regression.pdf"
#'
#' @return Plots to a pdf
#' 
#' @export
plotPairscan <- function(data.obj, pairscan.obj, phenotype = NULL, standardized = FALSE,
	show.marker.labels = FALSE, show.chr = TRUE, label.chr = TRUE, show.alleles = TRUE,
	allele.labels = NULL, pos.col = "brown", neg.col = "blue", 
	color.scheme = c("DO/CC", "other"), pdf.label = "Pairscan.Regression.pdf") {
  
  pairscan.results <- pairscan.obj$pairscan.results
  
  marker.pairs <- pairscan.results[[1]][[1]][,1:2]
  #get the markers used in the pair scan and sort them.
  markers <- unique(as.vector(marker.pairs))
  split.markers <- strsplit(markers, "_")
  markers.no.allele <- sapply(split.markers, function(x) x[1])
  markers.just.allele <- sapply(split.markers, function(x) x[2])
  
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
    
    allele.colors <- get.allele.colors(color.scheme, markers.just.allele)
    allele.labels <- allele.colors[,2]
    
    all.alleles <- unlist(lapply(strsplit(sorted.markers, "_"), function(x) x[2]))
    allele.cols <- allele.colors[match(all.alleles, all.alleles),3]
    
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