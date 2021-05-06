#' Identify linkage blocks
#' 
#' This function identifies linkage blocks among markers
#' using pairwise correlation between genotypes. 
#' The algorithm clusters adjacent, correlated markers
#' using the fastgreedy community detection algorithm from 
#' R/igraph
#' 
#' Csardi G, Nepusz T: The igraph software package for 
#' complex network research, InterJournal, Complex Systems
#' 1695. 2006. https://igraph.org
#' 
#' The correlation network can be optionally soft thresholded
#' to increase or decrease resolution. 
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' @param collapse_linked_markers A logical value. If TRUE markers are combined 
#' into linkage blocks based on correlation. If FALSE, each marker is treated as 
#' an independent observation.
#' @param threshold_power A soft threshold power. The marker
#' correlation matrix is raised to this power to increase or
#' decrease the number of linkage blocks detected. Increasing 
#' the power makes more linkage blocks, and decreasing the power
#' makes fewer linkage blocks. The default power is 1, which uses
#' the correlation matrix as is.
#' @param plot_blocks logical. If TRUE, the marker correlation 
#' matrices are plotted and the boundaries of the blocks are shown.
#' 
#' @return The data object is returned with a new list called
#' linkage_blocks_collapsed if collapse_linked_markers is TRUE 
#' and linkage_blocks_full if collapse_linked_markers is FALSE
#' Each element of the list is one linkage block and contains
#' a vector naming the markers in that block. Blocks are named
#' with a chromosome number and an index.
#' 
#' @keywords internal

linkage_blocks_network <- function(data_obj, geno_obj, collapse_linked_markers = TRUE, 
  threshold_power = 1, plot_blocks = TRUE){
  
  geno_names <- data_obj$geno_names
  marker_names <- geno_names[[3]]
  net_data <- data_obj$var_to_var_p_val
  
  
  if(length(net_data) == 0){
    stop("calc_p() must be run to calculate variant-to-variant influences.")
  }
  
  get_allele <- function(element){
    if(length(element) == 1){
      return(element[1])
    }else{
      return(element[2])
    }
  }
    
  get_marker_name <- function(element){
    return(element[1])
  }
  
  #get covariate information
  covar_info <- get_covar(data_obj)
  
  #find all the chromosomes that were used in the pairwise scan and sort them
  #with the refactoring geno_for_pairscan is no longer a reliable indicator 
  #of which markers were used in the pairscan.
  used_markers <- unique(as.vector(data_obj$var_to_var_p_val[,1:2]))
  all_marker_chr <- sapply(used_markers, function(x) get_marker_chr(data_obj, x))
  u_chr <- sort(as.numeric(unique(all_marker_chr)))

  if(u_chr[1] == 0){
    u_chr <- c(u_chr[-1], 0)
  }
  
  all_marker_names <- unlist(lapply(strsplit(used_markers, "_"), get_marker_name)) 
  marker_locale <- match(all_marker_names, geno_names[[3]])
  #========================================================================================
  # internal functions
  #========================================================================================
  # my_palette <- colorRampPalette(c("lightblue2", "green4"),space = "rgb")
  
  my_palette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  #if we are not collapsing the markers into blocks, 
  #or a chromosome only has one marker, or we are
  #adding the covariate chromosome, just add all 
  #individual markers to the list of blocks.
  add_ind_markers <- function(link_blocks, ch, chr_markers){
    chr_block_num <- 1
    if(is.null(link_blocks[[1]])){
      num_blocks = 1
    }else{
      num_blocks = length(link_blocks) + 1	
    }
    for(i in 1:length(chr_markers)){
      link_blocks[[num_blocks]] <- chr_markers[i]
      if(ch == 0){
        names(link_blocks)[num_blocks] <- chr_markers[i]
      }else{
        names(link_blocks)[num_blocks] <- paste("Chr", ch, "_", chr_block_num, sep = "")
      }
      num_blocks <- num_blocks + 1
      chr_block_num <- chr_block_num + 1
    }
    return(link_blocks)
  }
  
  
  #get the recombination data for the markers on a given chromosome
  get_chr_cor <- function(ordered_names, ordered_alleles){		
    marker_pos <- match(ordered_names, dimnames(geno_obj)[[3]])
    allele_pos <- match(ordered_alleles, dimnames(geno_obj)[[2]])
    chr_geno <- sapply(1:length(marker_pos), function(x) geno_obj[,allele_pos[x], marker_pos[x]])
    chr_cor <- cor(chr_geno, use = "complete.obs")
    return(chr_cor)
  }
  
  #========================================================================================
  # end internal functions
  #========================================================================================
  
  
  if(plot_blocks){pdf(paste("Recomb.Images.Genotype.Net.Thresh.", threshold_power, ".pdf", 
  sep = ""), width = 10, height = 5)}
  #go through each chromosome separately and find the linkage blocks on each chromosome
  link_blocks <- vector(mode = "list", length = 1)
  num_blocks <- 1
  for(ch in u_chr){
    chr_blocks = 1
    chr_markers <- used_markers[which(all_marker_chr == ch)]
    split_markers <- strsplit(chr_markers, "_")
    chr_alleles <- sapply(split_markers, function(x) x[2])
    chr_marker_names <- sapply(split_markers, function(x) x[1])
    
    # if(lookup.marker_position){
      # cat("looking up SNP positions...\n")
      # snp.info <- lapply(chr_marker_names, function(x) as.matrix(biomaRt::getBM(c("refsnp_id","allele","chr_name","chrom_start"), filters = "snp_filter", values = x, mart = snp.db)))
      # block.bp <- lapply(snp.info, function(x) as.numeric(x[,4]))
      # names(block.bp) <- chr_marker_names
      # no.info <- which(unlist(lapply(block.bp, length)) == 0)
      # block.bp[no.info] <- NA
      # block.bp <- unlist(block.bp)	
    # }else{
      #otherwise get positions from the data object
    block.bp <- get_marker_location(data_obj, chr_markers)
    # }
    
    if(!collapse_linked_markers || length(chr_markers) == 1 || ch == 0){
      link_blocks <- add_ind_markers(link_blocks, ch, chr_markers)
      # num_blocks <- num_blocks + length(link_blocks)
      num_blocks <- num_blocks + 1
    }else{
      marker_order <- order(block.bp)
      all_cor <- get_chr_cor(chr_marker_names[marker_order], chr_alleles[marker_order])
      diag(all_cor) <- 0
      thresh_mat <- abs(all_cor^threshold_power)
      net <- graph.adjacency(thresh_mat, mode = "undirected", weighted = TRUE)
      comm <- fastgreedy.community(net)$membership
      
      #In populations like the BXD, there is long-range LD that
      #complicates this blocking process. For now I will only 
      #block correlated markers that are also adjacent. This means
      #that two communities labeled 1 that are separated by other
      #communities, are actually two communities, and the second 
      #community 1 needs to have a new name.
      comm <- check_communities(comm)
      
      allele_table <- cbind(chr_markers, comm, chr_alleles)
      #sort by community and alleles so we don't break blocks
      #when alleles alternate back and forth. If we pick multiple
      #alleles for each marker, the markers will be in order, but
      #the alleles will cycle creating artificial block changes
      sorted_table <- sort_by_then_by(allele_table, sort_cols = c(2,3), col_type = c("n", "c"))
      
      chr_markers <- sorted_table[,1]
      comm <- as.numeric(sorted_table[,2])
      chr_alleles <- sorted_table[,3]
      
      allele_pairs <- consec_pairs(chr_alleles)
      allele_changes <- which(!apply(allele_pairs, 1, function(x) x[1] == x[2])) #find everywhere the community number 
      adj_comm <- consec_pairs(comm)
      cm_changes <- which(!apply(adj_comm, 1, function(x) x[1] == x[2])) #find everywhere the community number changes
      
      #each time the allele changes within a community, increment
      #that position and all the following positions
      allele_within_comm <- setdiff(allele_changes, cm_changes)
      if(length(allele_within_comm) > 0){
        for(ac in 1:length(allele_within_comm)){
          start_pos <- allele_within_comm[ac]+1
          end_pos <- length(comm)
          comm[start_pos:end_pos] <- comm[start_pos:end_pos] + 1
        }
      }
      
      #recalculate the community changes
      adj_comm <- consec_pairs(comm)
      cm_changes <- which(!apply(adj_comm, 1, function(x) x[1] == x[2])) #find everywhere the community number changes
      
      
      if(length(cm_changes) == 0){ #if there are no changes, put the whole chromosome into the block
        link_blocks[[num_blocks]] <- chr_markers
        names(link_blocks)[num_blocks] <- paste("Chr", ch, "_", chr_blocks, sep = "")
        num_blocks <- num_blocks + 1
      }else{ #otherwise, step through the communities and add each one as a block
        #making sure that different alleles are not grouped together
        
        #for each block on the chromosome
        for(cm in 1:(length(cm_changes)+1)){
          cm_locale <- which(comm == cm)
          marker_names <- chr_markers[cm_locale]
          
          link_blocks[[num_blocks]] <- marker_names
          names(link_blocks)[num_blocks] <- paste("Chr", ch, "_", chr_blocks, sep = "")
          num_blocks <- num_blocks + 1
          chr_blocks <- chr_blocks + 1
        } #end looping through communities
      } #end adding blocks based on communities
    } #end case for if we are not dealing with the covariate chromosome or not collapasing markers
    
    
    if(plot_blocks && ch != 0){
      
      layout(matrix(c(1,2), nrow = 1))
      image(1:dim(all_cor)[1], 1:dim(all_cor)[2], all_cor, main = paste("Marker Correlation Chr", ch), xlim = c(0,(dim(all_cor)[1]+1)), ylim = c(0,(dim(all_cor)[1]+1)), col = my_palette(50), axes = FALSE, ylab = "", xlab = "")
      
      block_names <- sapply(strsplit(names(link_blocks), "_"), function(x) x[1])
      chr_names <- sapply(strsplit(block_names, "Chr"), function(x) x[2])
      final_block_locale <- which(chr_names == ch)
      start_block = 0.5
      #outline each block
      for(b in 1:length(final_block_locale)){
        end_block <- start_block + length(link_blocks[[final_block_locale[b]]])
        segments(x0 = start_block, y0 = start_block, x1 = start_block, y1 = end_block, lwd = 3)
        segments(x0 = start_block, y0 = start_block, x1 = end_block, y1 = start_block, lwd = 3)
        segments(x0 = end_block, y0 = start_block, x1 = end_block, y1 = end_block, lwd = 3)
        segments(x0 = start_block, y0 = end_block, x1 = end_block, y1 = end_block, lwd = 3)
        start_block <- end_block
      } #end outlining blocks
      
      image(1:dim(thresh_mat)[1], 1:dim(thresh_mat)[2], thresh_mat, main = "Thresholded Correlations", xlim = c(0,(dim(thresh_mat)[1]+1)), ylim = c(0,(dim(thresh_mat)[1]+1)), col = my_palette(50), axes = FALSE, ylab = "", xlab = "")
      block_names <- sapply(strsplit(names(link_blocks), "_"), function(x) x[1])
      chr_names <- sapply(strsplit(block_names, "Chr"), function(x) x[2])
      final_block_locale <- which(chr_names == ch)
      start_block = 0.5
      #outline each block
      for(b in 1:length(final_block_locale)){
        end_block <- start_block + length(link_blocks[[final_block_locale[b]]])
        segments(x0 = start_block, y0 = start_block, x1 = start_block, y1 = end_block, lwd = 3)
        segments(x0 = start_block, y0 = start_block, x1 = end_block, y1 = start_block, lwd = 3)
        segments(x0 = end_block, y0 = start_block, x1 = end_block, y1 = end_block, lwd = 3)
        segments(x0 = start_block, y0 = end_block, x1 = end_block, y1 = end_block, lwd = 3)
        start_block <- end_block
      } #end outlining blocks		
      
    } #end plotting
    
    
  } #end looping through chromosomes	
  
  if(plot_blocks){dev.off()}
  
  if(collapse_linked_markers){
    data_obj$linkage_blocks_collapsed <- link_blocks
  }else{
    data_obj$linkage_blocks_full <- link_blocks
  }
  
  return(data_obj)
  
}