#' This script draws a head specifically over electrode locations
#' It shouldn't need to be scaled independently of the net layout
plotNetworkDO <- function(data.obj, marker.pairs = NULL, collapsed.net = TRUE, trait = NULL, 
                          phenotype.labels = NULL, ref.allele = "B", color.scheme = c("DO/CC", "other"),
                          main.lwd = 4, inter.lwd = 3, label.cex = 1.5, percent.bend = 15, chr.gap = 1, 
                          label.gap = 5, positive.col = "brown", negative.col = "blue", show.alleles = TRUE, 
                          allele.labels = NULL){
  
  if(collapsed.net){
    adj.mat <- data.obj$collapsed_net
    blocks <- data.obj$linkage_blocks_collapsed
  }else{
    adj.mat <- data.obj$full_net
    blocks <- data.obj$linkage_blocks_full
  }
  
  if(is.null(adj.mat)){
    stop("get.network() must be run before plotting the collapsed network.")
  }
  
  pos.col <- get.color(positive.col)[3]
  neg.col <- get.color(negative.col)[3]
  
  circle.dens = 0.0005
  center.x = 1; center.y = 1; radius = 2
  
  # chr.rel.length <- c(98.5, 103.9, 82.7, 88.6, 90.2, 79, 89.1, 76.2, 75.1, 77.9, 88, 63.9, 67.3, 66.4, 59, 57.8, 61.3, 59.4, 56.9)
  chr.rel.length <- c(1)
  
  all.chr <- data.obj$chromosome
  all.pos <- data.obj$marker_location
  all.pos[which(all.pos == 0)] <- 1 #we can't place markers as position 0. Change any 0s to 1.
  chr <- unique(all.chr)
  num.true.chr = length(chr)
  names(chr) <- rep("chr", num.true.chr)
  
  covar.info <- get.covar(data.obj)
  covar.names <- covar.info$covar.names
  if(length(covar.names) > 1){
    covar.pos1 <- segment.region(0, 0.9, length(covar.names))
    covar.pos2 <- segment.region(0.1, 1, length(covar.names))
  }else{
    covar.pos1 <- 0.5
    covar.pos2 <- 0.5	
  }
  
  if(length(covar.names) > 0){
    names(covar.names) <- rep("covar", length(covar.names))
    chr <- c(covar.names, chr)
  }
  
  if(length(chr.rel.length) != num.true.chr){
    # warning("The relative lengths of the chromosomes will not be plotted.")
    rel.length <- rep(1, num.true.chr)
  }else{
    rel.length <- chr.rel.length/max(chr.rel.length)
  }
  
  if(length(covar.names) > 0){
    rel.length <- c(rep(0.2, length(covar.names)), rel.length)
  }
  
  #============================================================================================
  # internal functions
  #============================================================================================
  
  get.block.coord <- function(radius.coord, start, pts.per.chr, block.rel.locale, idx, chr.blocks.locale, ch){
    coord.x <- radius.coord$x[start:(start+pts.per.chr[idx]-1)] #get the relative x and y coordinates for the block
    coord.y <- radius.coord$y[start:(start+pts.per.chr[idx]-1)]
    x.coord <- coord.x[round(length(coord.x)*as.numeric(block.rel.locale[2])):round(length(coord.x)*as.numeric(block.rel.locale[3]))]
    y.coord <- coord.y[round(length(coord.x)*as.numeric(block.rel.locale[2])):round(length(coord.y)*as.numeric(block.rel.locale[3]))]
    return(cbind(rep(names(chr.blocks.locale)[ch], length(x.coord)), x.coord, y.coord))
  }
  
  #assign a chromosome and relative position to each block
  get.chr.pos <- function(block){
    split.names <- strsplit(block, "_")
    just.locus <- unlist(lapply(split.names, function(x) x[1]))
    marker.locale <- which(data.obj$geno_names[[3]] %in% just.locus)
    if(length(marker.locale) == 0){ #if this marker isn't in the names, it's probably a covariate
      #look in the covariates
      marker.locale <- which(covar.names == block)
      return(c(0, covar.pos1[marker.locale], covar.pos2[marker.locale]))
    }
    chr <- unique(all.chr[marker.locale])
    if(length(chr) > 1){
      chr.char <- paste(chr, collapse = ", ")
      stop(paste("There is linkage between markers on chromosomes ", chr.char,". Please try a high r2.thresh.", sep = ""))
    }
    min.pos <- min(all.pos[marker.locale])
    max.pos <- max(all.pos[marker.locale])
    total.length <- max(all.pos[all.chr == chr])
    return(c(chr, min.pos/total.length, max.pos/total.length))
  }
  
  get.block.col <- function(block, allele.colors){
    markers <- blocks[[block]]
    #blocks are split by allele, so we only need to look at the first entry
    allele <- strsplit(markers, "_")[[1]][2]
    allele.locale <- which(allele.colors[,2] == allele)
    return(allele.colors[allele.locale,3])
  }
  #============================================================================================
  # end internal functions
  #============================================================================================
  
  #get coordinates for the concentric circles we will use
  chr.radius <- get.circle(radius, dens = circle.dens) #circle for chromosomes
  
  #these need to change if we are showing allele colors
  inner.bar.radius = get.circle(radius*0.98, dens = circle.dens) #circle for chr blocks if no alleles are shown
  
  #divide into chromosomes
  num.chr = length(chr)
  gap = round((length(chr.radius$x)*chr.gap)/100) #number of values to skip for gap between chromosomes
  
  
  label.gap <- round((length(chr.radius$x)*label.gap)/100)
  full.length <- length(chr.radius$x) - (gap*num.chr) - label.gap
  full.chr <- rel.length*(1/sum(rel.length))
  pts.per.chr <- floor(full.length*full.chr)
  
  if(!is.null(marker.pairs)){
    source.chr <- names(blocks)[sapply(marker.pairs[,1], function(x) grep(x, blocks))]
    target.chr <- names(blocks)[sapply(marker.pairs[,2], function(x) grep(x, blocks))]
    marker.pairs <- cbind(source.chr, target.chr)
  }else{
    #if no marker pairs are specified, take all of them
    just.int <- adj.mat[,1:nrow(adj.mat)]
    marker.pairs <- which(just.int != 0, arr.ind = TRUE)
    marker.pairs[,1] <- rownames(just.int)[as.numeric(marker.pairs[,1])]
    marker.pairs[,2] <- rownames(just.int)[as.numeric(marker.pairs[,2])]		
  }
  
  
  if(collapsed.net){
    blocks <- data.obj$linkage_blocks_collapsed
    #remove the blocks that were never tested
    not.tested <- setdiff(names(blocks), rownames(adj.mat))
    if(length(not.tested) > 0){
      not.tested.locale <- match(not.tested, names(blocks))
      blocks <- blocks[-not.tested.locale]
      data.obj$collapsed_net <- blocks
    }
  }else{
    blocks <- data.obj$linkage_blocks_full	
    #remove the blocks that were never tested
    not.tested <- setdiff(names(blocks), rownames(adj.mat))
    if(length(not.tested) > 0){
      not.tested.locale <- match(not.tested, names(blocks))
      blocks <- blocks[-not.tested.locale]
      data.obj$full_net <- blocks
    }
  }
  
  chr.pos <- t(sapply(blocks, get.chr.pos))
  colnames(chr.pos) <- c("chromosome", "min.position", "max.position")
  
  
  #and each phenotype
  if(is.null(trait)){
    pheno <- names(data.obj$max_var_to_pheno_influence)				
  }else{
    pheno <- trait
    trait.locale <- which(trait %in% names(data.obj$max_var_to_pheno_influence))
    if(length(trait.locale) < length(trait)){
      not.found <- which(!trait %in% names(data.obj$max_var_to_pheno_influence))
      message("I couldn't find the following traits:")
      cat(trait[not.found], sep = "\n")
      return()
    }
  }
  
  if(is.null(phenotype.labels)){
    pheno.names <- pheno
  }else{
    pheno.names <- phenotype.labels
    if(length(phenotype.labels) != length(pheno)){
      stop("I'm detecting the wrong number of phenotype labels.")
    }	
  }
  
  #get the circle coordinates for each trait
  #circles are generated from inner-most to
  #outer-most
  gap.rad = 0.05
  start.rad <- radius + gap.rad
  trait.circ <- get.concent.circ(pheno.names, start.rad = start.rad, gap.rad = gap.rad)
  new.start.rad <- start.rad + (gap.rad*length(pheno.names))
  
  #also get circles and colors for the different 
  #alleles if we are going to show them
  alleles <- unique(sapply(strsplit(colnames(data.obj$geno_for_pairscan), "_"), function(x) x[2]))
  if(show.alleles){
    allele.colors <- get.allele.colors(color.scheme, alleles)
  }
  
  
  label.radius = get.circle(new.start.rad+(gap.rad*3))
  
  #if we need to filter chr.pos and adj.mat to include
  #only the phenotypes we are including
  all.block.pheno <- c(rownames(chr.pos), pheno)
  adj.mat <- adj.mat[,colnames(adj.mat) %in% all.block.pheno, drop = FALSE]
  
  main.effect.mat <- adj.mat[,which(colnames(adj.mat) %in% pheno), drop = FALSE]
  
  chr.coord.table <- NULL #the coordinates of the chromosomes
  block.coord.table <- NULL #the coordinates of the blocks for plotting interaction polygons
  inner.bar.coord.table <- NULL #the coordinates of the blocks for plotting target bars
  
  plot.new()
  #give the right margin a bit more room to write covariate names
  plot.window(xlim = c(min(label.radius$x), max(label.radius$x)*1.25), ylim = c(min(label.radius$y), max(label.radius$y)))
  par(mar = c(2,2,2,0))
  plot.dim <- par("usr")
  
  #segment the y coordinates of the gap region
  #in the outermost circle to evenly space the
  #label sticks and gaps
  
  #add trait circles		
  plot.trait.circ(trait.circ, label.gap, plot.dim, main.lwd)
  
  
  start = label.gap
  for(i in 1:length(chr)){
    chr.x.coord <- chr.radius$x[start:(start+pts.per.chr[i]-1)]
    chr.y.coord <- chr.radius$y[start:(start+pts.per.chr[i]-1)]
    points(chr.x.coord, chr.y.coord, type = "l", lwd = main.lwd)
    chr.coord.table <- rbind(chr.coord.table, cbind(rep(chr[i], length(chr.x.coord)), chr.x.coord, chr.y.coord))
    
    if(names(chr)[i] == "covar"){
      text.adj = 0
    }else{
      text.adj = 0.5
    }
    
    text(mean(label.radius$x[start:(start+pts.per.chr[i]-1)]), mean(label.radius$y[start:(start+pts.per.chr[i]-1)]), chr[i], adj = text.adj, cex = label.cex)
    
    # get the number of blocks on this chromosome
    if(names(chr)[i] == "covar"){
      chr.blocks.locale <- which(rownames(chr.pos) == chr[i])
      names(chr.blocks.locale) <- chr[i]
    }else{
      chr.blocks.locale <- which(chr.pos[,1] == chr[i])
    }
    
    if(length(chr.blocks.locale) > 0){
      for(ch in 1:length(chr.blocks.locale)){
        main.effects <- main.effect.mat[chr.blocks.locale[ch],,drop=FALSE]
        block.rel.locale <- chr.pos[chr.blocks.locale[ch],,drop=FALSE]
        
        for(ph in 1:length(pheno)){
          
          if(main.effects[ph] != 0){ #if there are significant effects of this block, add them to the circle
            if(show.alleles && names(chr)[i] == "chr"){
              trait.col <- get.block.col(names(chr.blocks.locale)[ch], allele.colors = allele.colors)
            }else{
              if(main.effects[ph] < 0){trait.col = neg.col}else{trait.col = pos.col}
            }
            
            block.coord <- get.block.coord(radius.coord = trait.circ[[ph]], start, pts.per.chr, block.rel.locale, i, chr.blocks.locale, ch)
            if(nrow(block.coord) > 1){
              points(as.numeric(block.coord[,2]), as.numeric(block.coord[,3]), type = "l", lwd = main.lwd, col = trait.col)
            }else{
              points(as.numeric(block.coord[,2]), as.numeric(block.coord[,3]), type = "p", pch = 16, col = trait.col, cex = 0.7)	
            }
          }
          
          #collect positions of the blocks for polygons and inner target bars on slightly smaller circles
          block.coord <- get.block.coord(inner.bar.radius, start, pts.per.chr, block.rel.locale, i, chr.blocks.locale, ch)
          block.coord.table <- rbind(block.coord.table, block.coord)
          inner.bar.coord.table <- rbind(inner.bar.coord.table, block.coord)
          
        } #end looping through phenotypes
      } #end looping through blocks
    } #end case for if there are blocks on this chromosome
    start = start + pts.per.chr[i] + gap
  }
  
  
  #add the interactions
  just.m <- adj.mat[,-which(colnames(adj.mat) %in% pheno), drop = FALSE]
  
  if(!is.null(just.m)){
    new.mat <- matrix(0, nrow = nrow(just.m), ncol = ncol(just.m))
    colnames(new.mat) <- colnames(just.m)
    rownames(new.mat) <- rownames(just.m)
    source.ind <- match(marker.pairs[,1], rownames(new.mat))
    target.ind <- match(marker.pairs[,2], colnames(new.mat))
    for(i in 1:nrow(marker.pairs)){
      new.mat[source.ind[i], target.ind[i]] <- just.m[source.ind[i], target.ind[i]]
    }
    just.m <- new.mat
  }
  
  
  for(i in 1:nrow(just.m)){
    sig.locale <- which(just.m[i,] != 0)
    
    if(length(sig.locale) > 0){
      
      for(s in 1:length(sig.locale)){
        if(just.m[i,sig.locale[s]] > 0){edge.col <- pos.col}else{edge.col = neg.col}
        
        #find the start and stop positions
        start.block <- rownames(just.m)[i]
        start.block.coord <- block.coord.table[which(block.coord.table[,1] == start.block),,drop=FALSE]
        
        end.block <- colnames(just.m)[sig.locale[s]]
        end.block.coord <- block.coord.table[which(block.coord.table[,1] == end.block),,drop=FALSE]
        
        # #figure out the allele color(s) for the blocks
        # start.block.alleles <- strsplit(start.block, "_")[[1]]
        # end.block.alleles <- strsplit(end.block, "_")[[1]]
        
        # start.block.id <- match(start.block.alleles, alleles)
        # end.block.id <- match(end.block.alleles, alleles)
        
        #draw a polygon to connect the start block and stop position
        start.inter.x <- as.numeric(start.block.coord[,2])
        start.inter.y <- as.numeric(start.block.coord[,3])
        
        end.inter.x <- as.numeric(end.block.coord[,2])
        end.inter.y <- as.numeric(end.block.coord[,3])
        
        #sort the corners of the polygon according to their distance from the center of the circle
        start.mat <- matrix(c(start.inter.x[1], start.inter.y[1], start.inter.x[length(start.inter.x)], start.inter.y[length(start.inter.y)]), ncol = 2, byrow = TRUE)
        end.mat <- matrix(c(end.inter.x[1], end.inter.y[1], end.inter.x[length(end.inter.x)], end.inter.y[length(end.inter.y)]), ncol = 2, byrow = TRUE)
        rownames(start.mat) <- c("dist.start.min", "dist.start.max") 
        rownames(end.mat) <- c("dist.end.min", "dist.end.max") 
        
        #add a bar to indicate the block at the source end of the interaction
        start.bar.coord <- inner.bar.coord.table[which(inner.bar.coord.table[,1] == start.block),,drop=FALSE]
        
        points(as.numeric(start.bar.coord[,2]), as.numeric(start.bar.coord[,3]), type = "l", lwd = main.lwd, col = "darkgray")
        
        #add a bar to indicate the block at the target end of the interaction
        end.bar.coord <- inner.bar.coord.table[which(inner.bar.coord.table[,1] == end.block),,drop=FALSE]
        
        points(as.numeric(end.bar.coord[,2]), as.numeric(end.bar.coord[,3]), type = "l", lwd = main.lwd, col = "darkgray")
        
        
        #find the midpoint of the line between source and target
        mid.point.x <- mean(c(mean(start.inter.x), mean(end.inter.x)))
        mid.point.y <- mean(c(mean(start.inter.y), mean(end.inter.y)))
        
        
        #get the points on the line between the midpoint and the center of the circle
        pts.to.center <- get.line(mid.point.x, mid.point.y, center.x, center.y, dens = circle.dens)
        #find the point that makes this line the correct percentage long
        shifted.x <- pts.to.center$x[round(((percent.bend/100))*length(pts.to.center$x))]
        shifted.y <- pts.to.center$y[round(((percent.bend/100))*length(pts.to.center$y))]
        if(length(shifted.x) == 0){shifted.x = mid.point.x}
        if(length(shifted.y) == 0){shifted.y = mid.point.y}
        
        
        #make a spline curve based on the start and stop of the line plus this shifted midpoint
        inter.curve <- xspline(c(mean(start.inter.x), shifted.x, mean(end.inter.x)), y = c(mean(start.inter.y), shifted.y, mean(end.inter.y)), shape = 1, draw = FALSE)
        points(inter.curve$x, inter.curve$y, col = edge.col, type = "l", lwd = inter.lwd)
        
        #add an arrowhead pointed at the target
        arrow.rad <- atan2((mean(end.inter.y)-shifted.y), (mean(end.inter.x)-shifted.x))
        arrow.deg <- arrow.rad*(180/pi)
        shape::Arrowhead(x0 = mean(end.inter.x), y0 = mean(end.inter.y), arr.col = edge.col, arr.adj = 1, lcol = edge.col, angle = arrow.deg, arr.lwd = inter.lwd)
        
      }
    } #end case for when there are significicant interactions in this row
  }
  
  legend("bottomright", lty = 1, lwd = inter.lwd, col = c(pos.col, neg.col), legend = c("Positive", "Negative"))
  
  
}