#' This function uses the genotype correlations, but
#' finds blocks within the correlation matrix 
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param collapse.linked.markers A logical value. If TRUE markers are combined 
#' into linkage blocks based on correlation. If FALSE, each marker is treated as 
#' an independent observation.
#' @param threshold.power 
#' @param plot.blocks default = TRUE
#' @param lookup.marker.position default = TRUE
#' 
#' @param p.or.q A threshold indicating the maximum adjusted p value considered 
#' significant. If an fdr method has been used to correct for multiple testing, 
#' this value specifies the maximum q value considered significant.
#' @param min.std.effect 
#' @param standardize A logical value. If FALSE, the interaction terms are stored 
#' as marker to marker influence coefficients. If TRUE, the coefficients are 
#' standardized by their standard errors.
#' @param verbose default = FALSE
#' @param plot.linkage.blocks 
#' 
linkage.blocks.network <- function(data.obj, geno.obj, collapse.linked.markers = TRUE, 
threshold.power = 1, plot.blocks = TRUE, lookup.marker.position = FALSE){
  
  if(lookup.marker.position){
    require(biomaRt)
    if(data.obj$organism == "mouse"){		
      # lib <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
      lib <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "may2017.archive.ensembl.org")
      snp.db = biomaRt::useMart(biomart="ENSEMBL_MART_SNP", dataset="mmusculus_snp", host = "may2017.archive.ensembl.org")	
    }else{
      # lib <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")	
      lib <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "may2017.archive.ensembl.org")
      snp.db = biomaRt::useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp", host = "may2017.archive.ensembl.org")
    }
  }
  
  
  
  geno.names <- data.obj$geno_names
  marker.names <- geno.names[[3]]
  net.data <- data.obj$var_to_var_p_val
  
  
  if(length(net.data) == 0){
    stop("calc.p() must be run to calculate variant-to-variant influences.")
  }
  
  get.allele <- function(element){
    if(length(element) == 1){
      return(element[1])
    }else{
      return(element[2])
    }
  }
  
  get.chr <- function(element){
    if(length(element) == 1){
      return(0)
    }else{
      return(data.obj$chromosome[which(geno.names[[3]] == element[1])])
    }
  }
  
  get.marker.name <- function(element){
    return(element[1])
  }
  
  #get covariate information
  covar.info <- get.covar(data.obj)
  
  #find all the chromosomes that were used in the pairwise scan and sort them
  #with the refactoring geno_for_pairscan is no longer a reliable indicator 
  #of which markers were used in the pairscan.
  #used.markers <- colnames(data.obj$geno_for_pairscan) 
  used.markers <- unique(as.vector(data.obj$var_to_var_p_val[,1:2]))
  all.marker.chr <- unlist(lapply(strsplit(used.markers, "_"), get.chr))
  u_chr <- sort(as.numeric(unique(all.marker.chr)))

  if(u_chr[1] == 0){
    u_chr <- c(u_chr[-1], 0)
  }
  
  all.marker.names <- unlist(lapply(strsplit(used.markers, "_"), get.marker.name)) 
  marker.locale <- match(all.marker.names, geno.names[[3]])
  #========================================================================================
  # internal functions
  #========================================================================================
  # my.palette <- colorRampPalette(c("lightblue2", "green4"),space = "rgb")
  
  my.palette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  #if we are not collapsing the markers into blocks, 
  #or a chromosome only has one marker, or we are
  #adding the covariate chromosome, just add all 
  #individual markers to the list of blocks.
  add.ind.markers <- function(link.blocks, ch, chr.markers){
    chr.block.num <- 1
    if(is.null(link.blocks[[1]])){
      num.blocks = 1
    }else{
      num.blocks = length(link.blocks) + 1	
    }
    for(i in 1:length(chr.markers)){
      link.blocks[[num.blocks]] <- chr.markers[i]
      if(ch == 0){
        names(link.blocks)[num.blocks] <- chr.markers[i]
      }else{
        names(link.blocks)[num.blocks] <- paste("Chr", ch, "_", chr.block.num, sep = "")
      }
      num.blocks <- num.blocks + 1
      chr.block.num <- chr.block.num + 1
    }
    return(link.blocks)
  }
  
  
  #get the recombination data for the markers on a given chromosome
  get.chr.cor <- function(ordered.names, ordered.alleles){		
    marker.pos <- match(ordered.names, dimnames(geno.obj)[[3]])
    allele.pos <- match(ordered.alleles, dimnames(geno.obj)[[2]])
    chr.geno <- sapply(1:length(marker.pos), function(x) geno.obj[,allele.pos[x], marker.pos[x]])
    chr.cor <- cor(chr.geno, use = "complete.obs")
    return(chr.cor)
  }
  
  #========================================================================================
  # end internal functions
  #========================================================================================
  
  
  if(plot.blocks){pdf(paste("Recomb.Images.Genotype.Net.Thresh.", threshold.power, ".pdf", 
  sep = ""), width = 10, height = 5)}
  #go through each chromosome separately and find the linkage blocks on each chromosome
  link.blocks <- vector(mode = "list", length = 1)
  num.blocks <- 1
  for(ch in u_chr){
    chr.blocks = 1
    chr.markers <- used.markers[which(all.marker.chr == ch)]
    split.markers <- strsplit(chr.markers, "_")
    chr.alleles <- sapply(split.markers, function(x) x[2])
    chr.marker.names <- sapply(split.markers, function(x) x[1])
    
    if(lookup.marker.position){
      cat("looking up SNP positions...\n")
      snp.info <- lapply(chr.marker.names, function(x) as.matrix(biomaRt::getBM(c("refsnp_id","allele","chr_name","chrom_start"), filters = "snp_filter", values = x, mart = snp.db)))
      block.bp <- lapply(snp.info, function(x) as.numeric(x[,4]))
      names(block.bp) <- chr.marker.names
      no.info <- which(unlist(lapply(block.bp, length)) == 0)
      block.bp[no.info] <- NA
      block.bp <- unlist(block.bp)	
    }else{
      #otherwise get positions from the data object
      block.bp <- get.marker.location(data.obj, chr.markers)
    }
    
    if(!collapse.linked.markers || length(chr.markers) == 1 || ch == 0){
      link.blocks <- add.ind.markers(link.blocks, ch, chr.markers)
      # num.blocks <- num.blocks + length(link.blocks)
      num.blocks <- num.blocks + 1
    }else{
      marker.order <- order(block.bp)
      all.cor <- get.chr.cor(chr.marker.names[marker.order], chr.alleles[marker.order])
      diag(all.cor) <- 0
      thresh.mat <- abs(all.cor^threshold.power)
      net <- graph.adjacency(thresh.mat, mode = "undirected", weighted = TRUE)
      comm <- fastgreedy.community(net)$membership
      
      #In populations like the BXD, there is long-range LD that
      #complicates this blocking process. For now I will only 
      #block correlated markers that are also adjacent. This means
      #that two communities labeled 1 that are separated by other
      #communities, are actually two communities, and the second 
      #community 1 needs to have a new name.
      comm <- check.communities(comm)
      
      allele.table <- cbind(chr.markers, comm, chr.alleles)
      #sort by community and alleles so we don't break blocks
      #when alleles alternate back and forth. If we pick multiple
      #alleles for each marker, the markers will be in order, but
      #the alleles will cycle creating artificial block changes
      sorted.table <- sortByThenBy(allele.table, sort.cols = c(2,3), col.type = c("n", "c"))
      
      chr.markers <- sorted.table[,1]
      comm <- as.numeric(sorted.table[,2])
      chr.alleles <- sorted.table[,3]
      
      allele.pairs <- consec.pairs(chr.alleles)
      allele.changes <- which(!apply(allele.pairs, 1, function(x) x[1] == x[2])) #find everywhere the community number 
      adj.comm <- consec.pairs(comm)
      cm.changes <- which(!apply(adj.comm, 1, function(x) x[1] == x[2])) #find everywhere the community number changes
      
      #each time the allele changes within a community, increment
      #that position and all the following positions
      allele.within.comm <- setdiff(allele.changes, cm.changes)
      if(length(allele.within.comm) > 0){
        for(ac in 1:length(allele.within.comm)){
          start.pos <- allele.within.comm[ac]+1
          end.pos <- length(comm)
          comm[start.pos:end.pos] <- comm[start.pos:end.pos] + 1
        }
      }
      
      #recalculate the community changes
      adj.comm <- consec.pairs(comm)
      cm.changes <- which(!apply(adj.comm, 1, function(x) x[1] == x[2])) #find everywhere the community number changes
      
      
      if(length(cm.changes) == 0){ #if there are no changes, put the whole chromosome into the block
        link.blocks[[num.blocks]] <- chr.markers
        names(link.blocks)[num.blocks] <- paste("Chr", ch, "_", chr.blocks, sep = "")
        num.blocks <- num.blocks + 1
      }else{ #otherwise, step through the communities and add each one as a block
        #making sure that different alleles are not grouped together
        
        #for each block on the chromosome
        for(cm in 1:(length(cm.changes)+1)){
          cm.locale <- which(comm == cm)
          marker.names <- chr.markers[cm.locale]
          
          link.blocks[[num.blocks]] <- marker.names
          names(link.blocks)[num.blocks] <- paste("Chr", ch, "_", chr.blocks, sep = "")
          num.blocks <- num.blocks + 1
          chr.blocks <- chr.blocks + 1
        } #end looping through communities
      } #end adding blocks based on communities
    } #end case for if we are not dealing with the covariate chromosome or not collapasing markers
    
    
    if(plot.blocks && ch != 0){
      
      layout(matrix(c(1,2), nrow = 1))
      image(1:dim(all.cor)[1], 1:dim(all.cor)[2], all.cor, main = paste("Marker Correlation Chr", ch), xlim = c(0,(dim(all.cor)[1]+1)), ylim = c(0,(dim(all.cor)[1]+1)), col = my.palette(50), axes = FALSE, ylab = "", xlab = "")
      
      block.names <- sapply(strsplit(names(link.blocks), "_"), function(x) x[1])
      chr.names <- sapply(strsplit(block.names, "Chr"), function(x) x[2])
      final.block.locale <- which(chr.names == ch)
      start.block = 0.5
      #outline each block
      for(b in 1:length(final.block.locale)){
        end.block <- start.block + length(link.blocks[[final.block.locale[b]]])
        segments(x0 = start.block, y0 = start.block, x1 = start.block, y1 = end.block, lwd = 3)
        segments(x0 = start.block, y0 = start.block, x1 = end.block, y1 = start.block, lwd = 3)
        segments(x0 = end.block, y0 = start.block, x1 = end.block, y1 = end.block, lwd = 3)
        segments(x0 = start.block, y0 = end.block, x1 = end.block, y1 = end.block, lwd = 3)						
        start.block <- end.block
      } #end outlining blocks
      
      image(1:dim(thresh.mat)[1], 1:dim(thresh.mat)[2], thresh.mat, main = "Thresholded Correlations", xlim = c(0,(dim(thresh.mat)[1]+1)), ylim = c(0,(dim(thresh.mat)[1]+1)), col = my.palette(50), axes = FALSE, ylab = "", xlab = "")						
      block.names <- sapply(strsplit(names(link.blocks), "_"), function(x) x[1])
      chr.names <- sapply(strsplit(block.names, "Chr"), function(x) x[2])
      final.block.locale <- which(chr.names == ch)
      start.block = 0.5
      #outline each block
      for(b in 1:length(final.block.locale)){
        end.block <- start.block + length(link.blocks[[final.block.locale[b]]])
        segments(x0 = start.block, y0 = start.block, x1 = start.block, y1 = end.block, lwd = 3)
        segments(x0 = start.block, y0 = start.block, x1 = end.block, y1 = start.block, lwd = 3)
        segments(x0 = end.block, y0 = start.block, x1 = end.block, y1 = end.block, lwd = 3)
        segments(x0 = start.block, y0 = end.block, x1 = end.block, y1 = end.block, lwd = 3)						
        start.block <- end.block
      } #end outlining blocks		
      
    } #end plotting
    
    
  } #end looping through chromosomes	
  
  if(plot.blocks){dev.off()}
  
  if(collapse.linked.markers){
    data.obj$linkage_blocks_collapsed <- link.blocks
  }else{
    data.obj$linkage_blocks_full <- link.blocks
  }
  
  return(data.obj)
  
}