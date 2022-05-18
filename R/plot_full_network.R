#' Plot the final epistatic network in a traditional network view.
#' 
#'  This function plots the final results in a layout different to
#'  both \code{\link{plot_variant_influences}} and \code{\link{plot_network}}. 
#'  In this view, the network is plotted with a traditional network layout. 
#'  The genomic position information in \code{\link{plot_network}} is lost, but 
#'  in this view it is easier to see the structure of the overall network 
#'  in terms of hubs and peripheral nodes. In this view, each node is plotted 
#'  as a pie-chart, and the main effects of the node are indicated as 
#'  positive, negative, or not-significant (gray). Significant 
#'  interactions are shown arrows between 
#'  nodes and colored based on whether they are positive or negative interactions. 
#'  Colors for positive and negative main effects and interactions are specified
#' in the arguments. The function \code{\link{get_network}} must be run before plotting 
#'  the network.
#' 
#' For most networks, the default options will be fine, but there is a lot of room
#' for modification if changes are desired
#' 
#' @param data_obj A \code{\link{Cape}} object
#' @param p_or_q The maximum p value (or q value if FDR was used) for significant 
#' main effects and interactions.
#' @param collapsed_net A logical value indicating whether to show the network
#' condensed into linkage blocks (TRUE) or each individual marker (FALSE)
#' @param main A title for the plot
#' @param color_scheme either "DO/CC" or "other". "DO/CC" uses the official "DO/CC"
#' colors for the DO/CC alleles  
#' \url{http://www.csbio.unc.edu/CCstatus/index.py?run=AvailableLines.information}
#' "other" uses an unrelated color palette for multiple alleles.
#' @param pos_col The color to use for positive main effects and interactions
#' must be one of "green", "purple", "red", "orange", "blue", "brown", "yellow", "gray"
#' see \code{\link{get_color}}
#' @param neg_col The color to use for negative main effects and interactions
#' takes the same values as pos_col.
#' @param bg_col The color to be used in pie charts for non-significant main effects.
#' Takes the same values as pos_col
#' @param light_dark Indicates whether pos_col, neg_col, and bg_col should be selected
#' from light colors ("l"), dark colors ("d") or the full spectrum from light to dark ("f")
#' @param node_border_lwd The thickness of the lines around the pie charts
#' @param layout_matrix Users have the option of providing their own layout matrix for the
#' network. This should be a two column matrix indicating the x and y coordinates of each 
#' node in the network.
#' @param zoom Allows the user to zoom in and out on the image if the network is either 
#' running off the edges of the plot or too small in the middle of the plot.
#' @param xshift A constant by which to shift the x values of all nodes in the network.
#' @param yshift A constant by which to shift the y values of all nodes in the network.
#' @param node_radius The size of the pie chart for each node.
#' @param label_nodes A logical value indicating whether the nodes should be labeled.
#' Users may want to remove labels for large networks.
#' @param label_offset The amount by which to offset the node labels from the center of
#' the nodes.
#' @param label_cex The size of the node labels
#' @param legend_radius The size of the legend indicating which pie piece corresponds to which
#' traits.
#' @param legend_cex The size of the labels in the legend.
#' @param legend_position The position of the legend on the plot
#' @param arrow_offset The distance from the center of the node to the arrow head.
#' @param arrow_length The length of the head of the arrow
#' @param edge_lwd The thickness of the arrows showing the interactions.
#'
#' @return This function invisibly returns a list of length two. The first element
#' contains the igraph network object. The second contains the layout matrix for the
#' network. This can be passed in as an argument ("layout_matrix") which provides more
#' control to the user in the layout. Other network layouts from igraph can also be passed 
#' in here.
#'
#' @references Csardi G, Nepusz T: The igraph software package for complex network 
#' research, InterJournal, Complex Systems 1695. 2006. \url{https://igraph.org/}
#' 
#' @importFrom graphics arrows
#' @importFrom stats dist na.omit
#' 
#' @examples 
#' \dontrun{
#' plot_full_network(data_obj)
#' }
#' 
#' @export
 
plot_full_network <- function(data_obj, p_or_q = 0.05,  collapsed_net = TRUE, main = NULL, 
                         color_scheme = c("DO/CC", "other"), pos_col = "brown", neg_col = "blue", 
                         bg_col = "gray", light_dark = "f", node_border_lwd = 1, layout_matrix = NULL, 
                         zoom = 1, xshift = 0, yshift = 0, node_radius = 1, label_nodes = TRUE, 
                         label_offset = 0, label_cex = 1, legend_radius = 1, legend_cex = 1, 
                         legend_position = "topleft", arrow_offset = node_radius, arrow_length = 0.2, 
                         edge_lwd = 2){
  
  oldPar <- par(no.readonly = TRUE)
	on.exit(oldPar)

  pheno_tables <- data_obj$max_var_to_pheno_influence
  phenotypes <- names(pheno_tables)	
  
  
  plot_net_edges <- function(net, net_layout, lwd = 1, edge_col = "gray", arrow_offset = 1){
    
    edge_list <- get.edgelist(net, names = FALSE)
    
    if(length(col) == 1){
      edge_col = rep(edge_col, dim(edge_list)[1])
    }
    
    
    for(i in 1:length(edge_list[,1])){
      # cat(i, "\n")
      xy_start = c(net_layout[edge_list[i,1],1], net_layout[edge_list[i,1],2])
      xy_end <- c(net_layout[edge_list[i,2],1], net_layout[edge_list[i,2],2])
      alpha <- 0
      final_xy <- alpha*xy_start + (1-alpha)*(xy_end)
      dist_to_center <- dist(matrix(c(xy_end, final_xy), nrow = 2, byrow = TRUE), method = "euclidean")
      counter <- 0
      while(dist_to_center < arrow_offset && counter < 100){
        alpha <- alpha + 0.01
        final_xy <- alpha*xy_start + (1-alpha)*(xy_end)
        dist_to_center <- dist(matrix(c(xy_end, final_xy), nrow = 2, byrow = TRUE), method = "euclidean")	
        counter = counter + 1
        # print(counter)
      }
      # print(paste("\t", counter))
      # cat(round(xy_start, 2), " : ", round(final_xy, 2), "\n")
      arrows(x0 = xy_start[1], x1 = final_xy[1], y0 = xy_start[2], y1 = final_xy[2], lwd = lwd, col = edge_col[i], length = arrow_length)
    }
    
  }
  
  #This function plots nodes of a graph as polygons where each polygon
  #is subdivided into polygons each with a different color
  #This is for adding information to nodes about which phenotypes they
  #have main effects on
  plot_net_point <- function(x, y, node_radius = node_radius, 
  cols = c("green", "green", "red"), node_label = NULL, label_offset = 0, 
  label_cex = 1, border_col){
    
    draw_pie(x = x, y = y, radius = node_radius, cols = cols, add = TRUE, 
    node_border_lwd = node_border_lwd, border_col = border_col)
    
    
    if(!is.null(node_label)){
      text_x <- x+label_offset; text_y = y+label_offset
      text(x = text_x, y = text_y, labels = node_label, cex = label_cex)
    }
    
  }
  
  
  
  if(collapsed_net){
    blocks <- data_obj$linkage_blocks_collapsed
    net_data <- unlist(data_obj$collapsed_net)
    #convert this into an edge list to be compatible with the uncollapsed network
    sig_locale <- which(net_data != 0, arr.ind = TRUE)
    
    if(length(sig_locale) == 0){
      plot.new()
      plot.window(xlim = c(0,1), ylim = c(0,1))
      text(0.5, 0.5, "No Significant Interactions")
      invisible()
    }
    edge_vals <- net_data[which(net_data != 0)]
    sig_edges <- cbind(sig_locale, edge_vals)
    colnames(sig_edges) <- c("Source", "Target", "Effect")
    block_names <- matrix(NA, ncol = 2, nrow = dim(sig_edges)[1])
    block_names[,1] <- rownames(net_data)[sig_edges[,1]]
    block_names[,2] <- colnames(net_data)[sig_edges[,2]]
    sig_edges[,1:2] <- block_names
    
  }else{
    blocks = data_obj$linkage_blocks_full
    net_data <- data_obj$var_to_var_p_val
    
    var_sig_col <- which(colnames(data_obj$var_to_var_p_val) == "p_adjusted")
    sig_locale <- which(net_data[,var_sig_col] <= p_or_q)
    
    if(length(sig_locale) == 0){
      plot.new()
      plot.window(xlim = c(0,1), ylim = c(0,1))
      text(0.5, 0.5, "No Significant Interactions")
      invisible()
    }
    
    #get the edges for variant to variant interactions
    sig_edges <- matrix(net_data[sig_locale,1:2], nrow = length(sig_locale))
    effect_size <- (as.numeric(net_data[sig_locale,"Effect"])/as.numeric(net_data[sig_locale,"SE"]))
    sig_edges <- cbind(sig_edges, effect_size)
    colnames(sig_edges) <- c("Source", "Target", "Effect")
    
    #add any edges from the phenotype direct influences
    for(ph in 1:length(phenotypes)){
      pheno_stats <- pheno_tables[[ph]][,which(colnames(pheno_tables[[ph]]) != "t_stat")]
      p_adj_locale <- as.vector(na.omit(match(c("qval", "lfdr", "p_adjusted"), colnames(pheno_stats))))
      sig_ph_edges <- which(as.numeric(pheno_stats[,p_adj_locale]) <= p_or_q)
      if(length(sig_ph_edges) > 0){
        for(r in sig_ph_edges){
          m_coef <- as.numeric(pheno_stats[r,which(colnames(pheno_stats) == "coef")])
          m_se <- as.numeric(pheno_stats[r,which(colnames(pheno_stats) == "se")])
          table_row <- c(pheno_stats[r,1], names(pheno_tables)[ph], m_coef/m_se)
          sig_edges <- rbind(sig_edges, table_row)
        }
      }
    } #end adding phenotype edges
    
  }
  
  
  get_node_cols <- function(node_name, phenotypes, main_effects){
    node_locale <- which(main_effects[,"Source"] == node_name)
    no_effect_col <- bg_col
    node_cols <- rep(no_effect_col, length(phenotypes))
    
    node_effects <- main_effects[node_locale,c("Target", "Effect"),drop=FALSE]
    
    if(length(node_effects) > 0){
      for(i in 1:length(phenotypes)){
        pheno_locale <- which(node_effects[,"Target"] == phenotypes[i])
        if(length(pheno_locale) > 0){
          if(as.numeric(node_effects[pheno_locale,"Effect"]) > 0){node_cols[i] <- get_color(pos_col, light_dark)[2]}
          if(as.numeric(node_effects[pheno_locale,"Effect"]) < 0){node_cols[i] <- get_color(neg_col, light_dark)[2]}
        }
      }
    }
    return(node_cols)
  }
  
  
  main_effects_locale <- which(sig_edges[,2] %in% phenotypes)
  if(length(main_effects_locale) > 0){
    main_effects <- sig_edges[main_effects_locale,,drop=FALSE]
    interactions <- sig_edges[-main_effects_locale,,drop=FALSE]
  }else{
    main_effects <- NULL
    interactions <- sig_edges
  }
  
  
  if(length(interactions) == 0){
    plot.new()
    plot.window(xlim = c(0,1), ylim = c(0,1))
    text(0.5, 0.5, "No Significant Interactions")
    invisible()
  }else{	
    
    edgelist <- matrix(c(as.vector(interactions[,1]), as.vector(interactions[,2])), ncol = 2, byrow = FALSE)
        
    net <- graph.edgelist(edgelist)
    vertex_names <- V(net)$name
    
    if(!label_nodes){vertex_names <- NULL}
    
    edge_col <- rep(NA, ecount(net))
    edge_col[which(as.numeric(interactions[,"Effect"]) < 0)] <- get_color(neg_col, light_dark)[2]
    edge_col[which(as.numeric(interactions[,"Effect"]) > 0)] <- get_color(pos_col, light_dark)[2]
    
    if(is.character(layout_matrix)){
      if(layout_matrix == "manual"){
        tkp_id <- tkplot(net)
        readline(prompt = "Press return when ready:\n")
        layout_matrix <- tkplot.getcoords(tkp_id)
        tkplot.close(tkp_id)
      }else{
        layout_call <- call(layout_matrix, net)
        coord_matrix <- eval(layout_call)
      }			
    }else{ #otherwise layout_matrix is null or a matrix
      if(length(layout_matrix) == 0){
        coord_matrix <- layout.auto(net)
      }else{
        coord_matrix <- layout_matrix	
      }
    }
    #normalize the layout to be between -10 and 10
    norm_coord <- apply(coord_matrix, 2, function(x) x*10/max(abs(x)))
    
    #center the coordinate matrix, so we can spread out the nodes using zoom
    centered_coord <- apply(norm_coord, 2, function(x) x - mean(x))
    # minx <- min(centered_coord[,1]); maxx <- max(centered_coord[,1])
    # miny <- min(centered_coord[,2]); maxy <- max(centered_coord[,2])
    minx <- min(centered_coord); maxx <- max(centered_coord)
    miny <- min(centered_coord); maxy <- max(centered_coord)
    
    
    zoomed_coord <- centered_coord*zoom
    zoomed_coord[,1] <- zoomed_coord[,1] + xshift
    zoomed_coord[,2] <- zoomed_coord[,2] + yshift
    
    block_alleles <- sapply(V(net)$name, function(x) get_block_allele(data_obj, x, collapsed_net))
    allele_cols <- get_allele_colors(color_scheme, unique(block_alleles))
    block_cols <- allele_cols[match(block_alleles, allele_cols[,2]), 3]
    block_cols[which(is.na(block_cols))] <- "black"
    
    
    # par(mar = c(2,2,2,2))
    plot.new()
    plot.window(xlim = c(minx, maxx), ylim = c(miny, maxy))
    par(xpd = TRUE)
    plot_net_edges(net = net, net_layout = zoomed_coord, lwd = edge_lwd, 
    edge_col = edge_col, arrow_offset = arrow_offset)
    
    for(i in 1:length(V(net)$name)){
      plot_net_point(x = zoomed_coord[i,1], y = zoomed_coord[i,2], 
      node_radius = node_radius, 
      cols = get_node_cols(node_name = V(net)$name[i], phenotypes = phenotypes, main_effects = main_effects), 
      node_label = vertex_names[i], label_offset = label_offset, 
      label_cex = label_cex, border_col = block_cols[i])
    }
    
    if(legend_position == "topleft"){
      l_x <- minx; l_y <- maxy
    }
    if(legend_position == "topright"){
      l_x <- maxx; l_y <- maxy
    }
    if(legend_position == "bottomleft"){
      l_x <- minx; l_y <- miny
    }
    if(legend_position == "bottomright"){
      l_x <- maxx; l_y <- miny
    }
    
    
    draw_pie(x = l_x, y = l_y, radius = legend_radius, cols = rep(bg_col, length(phenotypes)), labels = phenotypes, label_cex = legend_cex, node_border_lwd = node_border_lwd, border_col = "black")
    
    if(!is.null(main)){
      text(x = mean(c(maxx, minx)), y = maxy, labels = main)
    }	
    
    par(xpd = FALSE)
    
    results <- list(net, zoomed_coord)
    names(results) <- c("net", "layout_matrix")
    invisible(results)
  }
  
  
}
