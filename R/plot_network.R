#' Plots cape results as a circular network
#'
#' This script plots cape results in a circular network.
#' The chromosomes are arranged in a circle. Main effects
#' are shown in concentric circles around the chromosomes,
#' with each trait in its own circle. Main effects can 
#' either be colored as negative or positive, or with parental
#' allele colors for multi-parent populations. 
#' 
#' Interaction effects are shown as arrows linking chromosomal
#' positions. They are colored based on whether they are positive
#' or negative.
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param marker_pairs a two-column matrix identifying which marker pairs should be plotted.
#' This is particularly useful if the network is very dense. The default value, NULL, plots
#' all marker pairs.
#' @param collapsed_net A logical value indicating whether to plot all individual SNPs
#' or linkage blocks calculated by \code{\link{linkage_blocks_network}}.
#' @param trait A character vector indicating which traits to plot. The default NULL
#' value plots all traits.
#' @param trait_labels A character vector indicating the names of the traits in case
#' the names from the data object are not great for plotting.
#' @param color_scheme A character value of either "DO/CC" or other indicating the 
#' color scheme of main effects. If "DO/CC" allele effects can be plotted with the
#' DO/CC colors.
#' @param main_lwd A numeric value indicating the line width for the main effect lines
#' @param inter_lwd A numeric value indicating the line width for the interaction lines
#' @param label_cex A numeric value indicating the size of the labels
#' @param percent_bend A numeric value indicating the amount that the arrows for the
#' interaction effects should be bent. A value of 0 will plot straight lines.
#' @param chr_gap A numeric value indicating the size of the gap plotted between chromosomes.
#' @param label_gap A numeric value indicating the size of the gap the chromosomes and their labels.
#' @param positive_col One of c("green", "purple", "red", "orange", "blue", "brown", "yellow", "gray")
#' indicating the color for positive interactions.
#' @param negative_col One of c("green", "purple", "red", "orange", "blue", "brown", "yellow", "gray")
#' indicating the color for negative interactions.
#' show_alleles A logical value indicating whether to color main effects by their allele values (TRUE)
#' or by whether they are positive or negative (FALSE)
#' @param show_alleles boolean, default is TRUE
#' 
#' @importFrom shape Arrowhead
#' @importFrom graphics xspline
#' 
#' @examples 
#' \dontrun{
#' plot_network(data_obj)
#' }
#' @export

plot_network <- function(data_obj, marker_pairs = NULL, collapsed_net = TRUE, 
                        trait = NULL, trait_labels = NULL, color_scheme = c("DO/CC", "other"), 
                        main_lwd = 4, inter_lwd = 3, label_cex = 1.5, percent_bend = 15, 
                        chr_gap = 1, label_gap = 5, positive_col = "brown", 
                        negative_col = "blue", show_alleles = TRUE){
  
    oldPar <- par(no.readonly = TRUE)
		on.exit(oldPar)

  if(collapsed_net){
    adj_mat <- data_obj$collapsed_net
    blocks <- data_obj$linkage_blocks_collapsed
  }else{
    adj_mat <- data_obj$full_net
    blocks <- data_obj$linkage_blocks_full
  }
  
  if(is.null(adj_mat)){
    stop("get_network() must be run before plotting the collapsed network.")
  }
  
  pos_col <- get_color(positive_col)[3]
  neg_col <- get_color(negative_col)[3]
  
  circle_dens = 0.0005
  center_x = 1; center_y = 1; radius = 2
  
  # chr_rel_length <- c(98.5, 103.9, 82.7, 88.6, 90.2, 79, 89.1, 76.2, 75.1, 77.9, 88, 63.9, 67.3, 66.4, 59, 57.8, 61.3, 59.4, 56.9)
  chr_rel_length <- c(1)
  
  all_chr <- data_obj$chromosome
  all_pos <- data_obj$marker_location
  all_pos[which(all_pos == 0)] <- 1 #we can't place markers as position 0. Change any 0s to 1.
  chr <- unique(all_chr)
  num_true_chr = length(chr)
  names(chr) <- rep("chr", num_true_chr)
  
  covar_info <- get_covar(data_obj)
  covar_names <- covar_info$covar_names
  if(length(covar_names) > 1){
    covar_pos1 <- segment_region(0, 0.9, length(covar_names))
    covar_pos2 <- segment_region(0.1, 1, length(covar_names))
  }else{
    covar_pos1 <- 0.5
    covar_pos2 <- 0.5	
  }
  
  if(length(covar_names) > 0){
    names(covar_names) <- rep("covar", length(covar_names))
    chr <- c(covar_names, chr)
  }
  
  if(length(chr_rel_length) != num_true_chr){
    # warning("The relative lengths of the chromosomes will not be plotted.")
    rel_length <- rep(1, num_true_chr)
  }else{
    rel_length <- chr_rel_length/max(chr_rel_length)
  }
  
  if(length(covar_names) > 0){
    rel_length <- c(rep(0.2, length(covar_names)), rel_length)
  }
  
  #============================================================================================
  # internal functions
  #============================================================================================
  
  get_block_coord <- function(radius_coord, start, pts_per_chr, block_rel_locale, idx, chr_blocks_locale, ch){
    coord_x <- radius_coord$x[start:(start+pts_per_chr[idx]-1)] #get the relative x and y coordinates for the block
    coord_y <- radius_coord$y[start:(start+pts_per_chr[idx]-1)]
    x_coord <- coord_x[round(length(coord_x)*as.numeric(block_rel_locale[2])):round(length(coord_x)*as.numeric(block_rel_locale[3]))]
    y_coord <- coord_y[round(length(coord_x)*as.numeric(block_rel_locale[2])):round(length(coord_y)*as.numeric(block_rel_locale[3]))]
    return(cbind(rep(names(chr_blocks_locale)[ch], length(x_coord)), x_coord, y_coord))
  }
  
  #assign a chromosome and relative position to each block
  get_chr_pos <- function(block){
    split_names <- strsplit(block, "_")
    just_locus <- unlist(lapply(split_names, function(x) x[1]))
    marker_locale <- which(data_obj$geno_names[[3]] %in% just_locus)
    if(length(marker_locale) == 0){ #if this marker isn't in the names, it's probably a covariate
      #look in the covariates
      marker_locale <- which(covar_names == block)
      return(c(0, covar_pos1[marker_locale], covar_pos2[marker_locale]))
    }
    chr <- unique(all_chr[marker_locale])
    if(length(chr) > 1){
      chr_char <- paste(chr, collapse = ", ")
      stop(paste("There is linkage between markers on chromosomes ", chr_char,". Please try a high r2.thresh.", sep = ""))
    }
    min_pos <- min(as.numeric(all_pos[marker_locale]))
    max_pos <- max(as.numeric(all_pos[marker_locale]))
    total_length <- max(as.numeric(all_pos[all_chr == chr]), na.rm =TRUE)
    return(c(chr, min_pos/total_length, max_pos/total_length))
  }
  
  get_block_col <- function(block, allele_colors){
    markers <- blocks[[block]]
    #blocks are split by allele, so we only need to look at the first entry
    allele <- strsplit(markers, "_")[[1]][2]
    allele_locale <- which(allele_colors[,2] == allele)
    return(allele_colors[allele_locale,3])
  }
  #============================================================================================
  # end internal functions
  #============================================================================================
  
  #get coordinates for the concentric circles we will use
  chr_radius <- get_circle(radius, dens = circle_dens) #circle for chromosomes
  
  #these need to change if we are showing allele colors
  inner_bar_radius = get_circle(radius*0.98, dens = circle_dens) #circle for chr blocks if no alleles are shown
  
  #divide into chromosomes
  num_chr = length(chr)
  gap = round((length(chr_radius$x)*chr_gap)/100) #number of values to skip for gap between chromosomes
  
  
  label_gap <- round((length(chr_radius$x)*label_gap)/100)
  full_length <- length(chr_radius$x) - (gap*num_chr) - label_gap
  full_chr <- rel_length*(1/sum(rel_length))
  pts_per_chr <- floor(full_length*full_chr)
  
  if(!is.null(marker_pairs)){
    source_chr <- names(blocks)[sapply(marker_pairs[,1], function(x) grep(x, blocks))]
    target_chr <- names(blocks)[sapply(marker_pairs[,2], function(x) grep(x, blocks))]
    marker_pairs <- cbind(source_chr, target_chr)
  }else{
    #if no marker pairs are specified, take all of them
    just_int <- adj_mat[,1:nrow(adj_mat)]
    marker_pairs <- which(just_int != 0, arr.ind = TRUE)
    marker_pairs[,1] <- rownames(just_int)[as.numeric(marker_pairs[,1])]
    marker_pairs[,2] <- rownames(just_int)[as.numeric(marker_pairs[,2])]		
  }
  
  if(collapsed_net){
    blocks <- data_obj$linkage_blocks_collapsed
    #remove the blocks that were never tested
    not_tested <- setdiff(names(blocks), rownames(adj_mat))
    if(length(not_tested) > 0){
      not_tested_locale <- match(not_tested, names(blocks))
      blocks <- blocks[-not_tested_locale]
      # TODO which line below is correct? The commented out line causes a bug in plot_full_network.R
      # TODO if the collapsed_net line is correct, then the setdiff on line 140 should be reversed
      # TODO and the collapsed_net object built again like it is in get_network line 152 .. 172
      #data_obj$collapsed_net <- blocks
      data_obj$linkage_blocks_collapsed <- blocks
    }
  }else{
    blocks <- data_obj$linkage_blocks_full	
    #remove the blocks that were never tested
    not_tested <- setdiff(names(blocks), rownames(adj_mat))
    if(length(not_tested) > 0){
      not_tested_locale <- match(not_tested, names(blocks))
      blocks <- blocks[-not_tested_locale]
      # TODO same as above...
      #data_obj$full_net <- blocks
      data_obj$linkage_blocks_full <- blocks
    }
  }
  
  chr_pos <- t(sapply(blocks, get_chr_pos))
  colnames(chr_pos) <- c("chromosome", "min_position", "max_position")
  
  
  #and each phenotype
  if(is.null(trait)){
    pheno <- names(data_obj$max_var_to_pheno_influence)				
  }else{
    pheno <- trait
    trait_locale <- which(trait %in% names(data_obj$max_var_to_pheno_influence))
    if(length(trait_locale) < length(trait)){
      not_found <- which(!trait %in% names(data_obj$max_var_to_pheno_influence))
      warning("I couldn't find the following traits:", paste(trait[not_found], collapse = ", "))
      stop()
    }
  }
  
  if(is.null(trait_labels)){
    pheno_names <- pheno
  }else{
    pheno_names <- trait_labels
    if(length(trait_labels) != length(pheno)){
      stop("I'm detecting the wrong number of phenotype labels.")
    }	
  }
  
  #get the circle coordinates for each trait
  #circles are generated from inner-most to
  #outer-most
  gap_rad = 0.05
  start_rad <- radius + gap_rad
  trait_circ <- get_concent_circ(pheno_names, start_rad = start_rad, gap_rad = gap_rad)
  new_start_rad <- start_rad + (gap_rad*length(pheno_names))
  
  #also get circles and colors for the different 
  #alleles if we are going to show them
  alleles <- unique(sapply(strsplit(rownames(adj_mat), "_"), function(x) x[2]))
  if(show_alleles){
    allele_colors <- get_allele_colors(color_scheme, alleles)
  }
  
  
  label_radius = get_circle(new_start_rad+(gap_rad*3))
  
  #if we need to filter chr_pos and adj_mat to include
  #only the phenotypes we are including
  all_block_pheno <- c(rownames(chr_pos), pheno)
  adj_mat <- adj_mat[,colnames(adj_mat) %in% all_block_pheno, drop = FALSE]
  
  main_effect_mat <- adj_mat[,which(colnames(adj_mat) %in% pheno), drop = FALSE]
  
  chr_coord_table <- NULL #the coordinates of the chromosomes
  block_coord_table <- NULL #the coordinates of the blocks for plotting interaction polygons
  inner_bar_coord_table <- NULL #the coordinates of the blocks for plotting target bars
  
  plot.new()
  #give the right margin a bit more room to write covariate names
  plot.window(xlim = c(min(label_radius$x), max(label_radius$x)*1.25), ylim = c(min(label_radius$y), max(label_radius$y)))
  par(mar = c(2,2,2,0))
  plot_dim <- par("usr")
  
  #segment the y coordinates of the gap region
  #in the outermost circle to evenly space the
  #label sticks and gaps
  
  #add trait circles		
  plot_trait_circ(trait_circ, label_gap, plot_dim, main_lwd)
  
  
  start = label_gap
  for(i in 1:length(chr)){
    chr_x_coord <- chr_radius$x[start:(start+pts_per_chr[i]-1)]
    chr_y_coord <- chr_radius$y[start:(start+pts_per_chr[i]-1)]
    points(chr_x_coord, chr_y_coord, type = "l", lwd = main_lwd)
    chr_coord_table <- rbind(chr_coord_table, cbind(rep(chr[i], length(chr_x_coord)), chr_x_coord, chr_y_coord))
    
    if(names(chr)[i] == "covar"){
      text_adj = 0
    }else{
      text_adj = 0.5
    }
    
    text(mean(label_radius$x[start:(start+pts_per_chr[i]-1)]), mean(label_radius$y[start:(start+pts_per_chr[i]-1)]), chr[i], adj = text_adj, cex = label_cex)
    
    # get the number of blocks on this chromosome
    if(names(chr)[i] == "covar"){
      chr_blocks_locale <- which(rownames(chr_pos) == chr[i])
      names(chr_blocks_locale) <- chr[i]
    }else{
      chr_blocks_locale <- which(chr_pos[,1] == chr[i])
    }
    
    if(length(chr_blocks_locale) > 0){
      for(ch in 1:length(chr_blocks_locale)){
        main_effects <- main_effect_mat[chr_blocks_locale[ch],,drop=FALSE]
        block_rel_locale <- chr_pos[chr_blocks_locale[ch],,drop=FALSE]
        
        for(ph in 1:length(pheno)){
          
          if(main_effects[ph] != 0){ #if there are significant effects of this block, add them to the circle
            if(show_alleles && names(chr)[i] == "chr"){
              trait_col <- get_block_col(names(chr_blocks_locale)[ch], allele_colors = allele_colors)
            }else{
              if(main_effects[ph] < 0){trait_col = neg_col}else{trait_col = pos_col}
            }
            
            block_coord <- get_block_coord(radius_coord = trait_circ[[ph]], start, pts_per_chr, block_rel_locale, i, chr_blocks_locale, ch)
            if(nrow(block_coord) > 1){
              points(as.numeric(block_coord[,2]), as.numeric(block_coord[,3]), type = "l", lwd = main_lwd, col = trait_col)
            }else{
              points(as.numeric(block_coord[,2]), as.numeric(block_coord[,3]), type = "p", pch = 16, col = trait_col, cex = 0.7)	
            }
          }
          
          #collect positions of the blocks for polygons and inner target bars on slightly smaller circles
          block_coord <- get_block_coord(inner_bar_radius, start, pts_per_chr, block_rel_locale, i, chr_blocks_locale, ch)
          block_coord_table <- rbind(block_coord_table, block_coord)
          inner_bar_coord_table <- rbind(inner_bar_coord_table, block_coord)
          
        } #end looping through phenotypes
      } #end looping through blocks
    } #end case for if there are blocks on this chromosome
    start = start + pts_per_chr[i] + gap
  }
  
  
  #add the interactions
  just_m <- adj_mat[,-which(colnames(adj_mat) %in% pheno), drop = FALSE]
  
  if(!is.null(just_m)){
    new_mat <- matrix(0, nrow = nrow(just_m), ncol = ncol(just_m))
    colnames(new_mat) <- colnames(just_m)
    rownames(new_mat) <- rownames(just_m)
    source_ind <- match(marker_pairs[,1], rownames(new_mat))
    target_ind <- match(marker_pairs[,2], colnames(new_mat))
    for(i in 1:nrow(marker_pairs)){
      new_mat[source_ind[i], target_ind[i]] <- just_m[source_ind[i], target_ind[i]]
    }
    just_m <- new_mat
  }
  
  
  for(i in 1:nrow(just_m)){
    sig_locale <- which(just_m[i,] != 0)
    
    if(length(sig_locale) > 0){
      
      for(s in 1:length(sig_locale)){
        if(just_m[i,sig_locale[s]] > 0){edge_col <- pos_col}else{edge_col = neg_col}
        
        #find the start and stop positions
        start_block <- rownames(just_m)[i]
        start_block_coord <- block_coord_table[which(block_coord_table[,1] == start_block),,drop=FALSE]
        
        end_block <- colnames(just_m)[sig_locale[s]]
        end_block_coord <- block_coord_table[which(block_coord_table[,1] == end_block),,drop=FALSE]
        
        # #figure out the allele color(s) for the blocks
        # start_block_alleles <- strsplit(start_block, "_")[[1]]
        # end_block_alleles <- strsplit(end_block, "_")[[1]]
        
        # start_block.id <- match(start_block_alleles, alleles)
        # end_block.id <- match(end_block_alleles, alleles)
        
        #draw a polygon to connect the start block and stop position
        start_inter_x <- as.numeric(start_block_coord[,2])
        start_inter_y <- as.numeric(start_block_coord[,3])
        
        end_inter_x <- as.numeric(end_block_coord[,2])
        end_inter_y <- as.numeric(end_block_coord[,3])
        
        #sort the corners of the polygon according to their distance from the center of the circle
        start_mat <- matrix(c(start_inter_x[1], start_inter_y[1], start_inter_x[length(start_inter_x)], start_inter_y[length(start_inter_y)]), ncol = 2, byrow = TRUE)
        end_mat <- matrix(c(end_inter_x[1], end_inter_y[1], end_inter_x[length(end_inter_x)], end_inter_y[length(end_inter_y)]), ncol = 2, byrow = TRUE)
        rownames(start_mat) <- c("dist.start.min", "dist.start.max") 
        rownames(end_mat) <- c("dist.end.min", "dist.end.max") 
        
        #add a bar to indicate the block at the source end of the interaction
        start_bar_coord <- inner_bar_coord_table[which(inner_bar_coord_table[,1] == start_block),,drop=FALSE]
        
        points(as.numeric(start_bar_coord[,2]), as.numeric(start_bar_coord[,3]), type = "l", lwd = main_lwd, col = "darkgray")
        
        #add a bar to indicate the block at the target end of the interaction
        end_bar_coord <- inner_bar_coord_table[which(inner_bar_coord_table[,1] == end_block),,drop=FALSE]
        
        points(as.numeric(end_bar_coord[,2]), as.numeric(end_bar_coord[,3]), type = "l", lwd = main_lwd, col = "darkgray")
        
        
        #find the midpoint of the line between source and target
        mid_point_x <- mean(c(mean(start_inter_x), mean(end_inter_x)))
        mid_point_y <- mean(c(mean(start_inter_y), mean(end_inter_y)))
        
        
        #get the points on the line between the midpoint and the center of the circle
        pts_to_center <- get_line(mid_point_x, mid_point_y, center_x, center_y, dens = circle_dens)
        #find the point that makes this line the correct percentage long
        shifted_x <- pts_to_center$x[round(((percent_bend/100))*length(pts_to_center$x))]
        shifted_y <- pts_to_center$y[round(((percent_bend/100))*length(pts_to_center$y))]
        if(length(shifted_x) == 0){shifted_x = mid_point_x}
        if(length(shifted_y) == 0){shifted_y = mid_point_y}
        
        
        #make a spline curve based on the start and stop of the line plus this shifted midpoint
        inter_curve <- xspline(c(mean(start_inter_x), shifted_x, mean(end_inter_x)), y = c(mean(start_inter_y), shifted_y, mean(end_inter_y)), shape = 1, draw = FALSE)
        points(inter_curve$x, inter_curve$y, col = edge_col, type = "l", lwd = inter_lwd)
        
        #add an arrowhead pointed at the target
        arrow_rad <- atan2((mean(end_inter_y)-shifted_y), (mean(end_inter_x)-shifted_x))
        arrow_deg <- arrow_rad*(180/pi)
        shape::Arrowhead(x0 = mean(end_inter_x), y0 = mean(end_inter_y), arr.col = edge_col, arr.adj = 1, lcol = edge_col, angle = arrow_deg, arr.lwd = inter_lwd)
        
      }
    } #end case for when there are significicant interactions in this row
  }
  
  legend("bottomright", lty = 1, lwd = inter_lwd, col = c(pos_col, neg_col), legend = c("Positive", "Negative"))
  
  
}