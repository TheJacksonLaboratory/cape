#' Plot the result of the pairwise scan
#'
#' This function plots the results of the pairwise scan.
#' It plots a matrix of the the interactions between all 
#' pairs of markers.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param pairscan_obj a pairscan object from \code{\link{pairscan}}
#' @param phenotype The names of the phenotypes to be plotted. If NULL,
#' all phenotypes are plotted.
#' @param standardized If TRUE, the standardized effects are plotted.
#' IF FALSE, the effect sizes are plotted.
#' @param show_marker_labels If TRUE marker labels are plotted along the
#' axes. If FALSE, they are omitted.
#' @param show_chr If TRUE, the chromosome boundaries are shown
#' @param label_chr If TRUE, the chromosomes are labeled
#' @param show_alleles If TRUE, the allele of each marker is indicated by color.
#' @param allele_labels Labels for the alleles if other than those stored in the
#' data object.
#' @param pos_col The color to use for positive main effects and interactions
#' must be one of "green", "purple", "red", "orange", "blue", "brown", "yellow", "gray"
#' see \code{\link{get_color}}
#' @param neg_col The color to use for negative main effects and interactions
#' takes the same values as pos_col.
#' @param color_scheme either "DO/CC" or "other". "DO/CC" uses the official "DO/CC"
#' colors for the DO/CC alleles  
#' \url{http://www.csbio.unc.edu/CCstatus/index.py?run=AvailableLines.information}
#' "other" uses an unrelated color palette for multiple alleles.
#' @param pdf_label Label for the resulting file. Defaults to "Pairscan.Regression.pdf"
#' if plotting to pdf, "Pairscan.Regression.jpg" otherwise.
#'
#' @return Plots to a pdf
#' 
#' @examples 
#' \dontrun{
#' plot_pairscan(data_obj, pairscan_obj)
#' }
#' 
#' @importFrom grDevices jpeg
#' 
#' @export
plot_pairscan <- function(data_obj, pairscan_obj, phenotype = NULL, standardized = FALSE,
	show_marker_labels = FALSE, show_chr = TRUE, label_chr = TRUE, show_alleles = TRUE,
	allele_labels = NULL, pos_col = "brown", neg_col = "blue", 
	color_scheme = c("DO/CC", "other"), pdf_label = "Pairscan.Regression.pdf") {
  
  pairscan_results <- pairscan_obj$pairscan_results
  
  marker_pairs <- pairscan_results[[1]][[1]][,1:2]
  #get the markers used in the pair scan and sort them.
  markers <- unique(as.vector(marker_pairs))
  split_markers <- strsplit(markers, "_")
  markers_no_allele <- sapply(split_markers, function(x) x[1])
  markers_just_allele <- sapply(split_markers, function(x) x[2])
  
  sorted_markers <- markers[order(as.numeric(get_marker_num(data_obj, markers_no_allele)))]
  
  #get coordinates of the chromosome boundaries
  if(show_chr){
    chromosomes <- get_marker_chr(data_obj, markers_no_allele)
    u_chr <- unique(chromosomes[!is.na(chromosomes)])
    chr_boundaries <- apply(matrix(u_chr, ncol = 1), 1, function(x) max(which(chromosomes == x))) + 0.5
    chr_boundaries <- c(0, chr_boundaries)
    if(label_chr){
      chr_names <- unique(chromosomes)
    }else{
      chr_names <- NULL	
    }
  }else{
    chr_boundaries <- NULL
    chr_names <- NULL
  }
  
  
  if(show_alleles){
    
    allele_colors <- get_allele_colors(color_scheme, markers_just_allele)
    allele_labels <- allele_colors[,2]
    
    all_alleles <- unlist(lapply(strsplit(sorted_markers, "_"), function(x) x[2]))
    allele_cols <- allele_colors[match(all_alleles, all_alleles),3]
    
  }else{
    allele_cols <- NULL
  }
  
  pairscan_result <- pairscan_obj$pairscan_results
  
  if(is.null(pairscan_result)){
    stop("pairscan() must be run before plot_pairscan()")
  }
  
  if(is.null(phenotype)){
    phenotype <- names(pairscan_result)
  }
  
  pheno_num <- which(names(pairscan_result) %in% phenotype)
  
  if(length(pheno_num) < length(phenotype)){
    not_found <- which(!(phenotype %in% names(pairscan_result)))
    warning("I couldn't find the following phenotypes: ", paste(phenotype[not_found], collapse = ", "))
    stop()
  }
  
  #collect the results, so we can put them on the same scale
  all_results_mats <- list()
  min_x <- 0
  max_x <- 0
  #for each phenotype scanned
  for(p in pheno_num){
    #build a results matrix
    results_mat <- matrix(0, length(markers), length(markers))
    colnames(results_mat) <- rownames(results_mat) <- sorted_markers
    
    if(standardized){
      pair_int <- as.numeric(pairscan_result[[p]][[1]][, "marker1:marker2"])/as.numeric(pairscan_result[[p]][[2]][,"marker1:marker2"])
    }else{
      pair_int <- as.numeric(pairscan_result[[p]][[1]][,"marker1:marker2"])
    }
    
    #create symetric matrices with the values
    net <- graph.edgelist(pairscan_result[[p]][[1]][,1:2])
    E(net)$weight <- pair_int
    results_mat <- as.matrix(as_adjacency_matrix(net, attr = "weight"))
    results_mat <- results_mat + t(results_mat)
    
    diag(results_mat) <- NA
    all_results_mats[[p]] <- results_mat
    min_x <- max(abs(c(max_x, results_mat)), na.rm = TRUE)*-1
    max_x <- max(abs(c(max_x, results_mat)), na.rm = TRUE)
    
  }	
  
  if (endsWith(pdf_label, "pdf")) {
    pdf(pdf_label, width = 7, height = 6)
  } else if (endsWith(pdf_label, "jpg")) {
    jpeg(pdf_label, quality = 100)
  }
  
  for(p in 1:length(pheno_num)){
    my_image_plot(x = all_results_mats[[pheno_num[p]]], xlab = "marker1", ylab = "marker2", main = phenotype[p], min_x = min_x, max_x = max_x, show_labels = show_marker_labels, chromosome_coordinates = chr_boundaries, chr_names = chr_names, allele_cols = allele_cols, pos_col = pos_col, neg_col = neg_col)
  }
  dev.off()
  
}