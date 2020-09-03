#' Plot the final epistatic network in a traditional network view.
#' 
#'  This function plots the final results in a layout different to
#'  both \code{\link{plotVariantInfluences}} and \code{\link{plotNetwork}}. 
#'  In this view, the network is plotted with a traditional network layout. 
#'  The genomic position information in \code{\link{plotNetwork}} is lost, but 
#'  in this view it is easier to see the structure of the overall network 
#'  in terms of hubs and peripheral nodes. In this view, each node is plotted 
#'  as a pie-chart, and the main effects of the node are indicated as 
#'  positive, negative, or not-significant (gray). Significant 
#'  interactions are shown arrows between 
#'  nodes and colored based on whether they are positive or negative interactions. 
#'  Colors for positive and negative main effects and interactions are specified
#' in the arguments. The function \code{\link{get.network}} must be run before plotting 
#'  the network.
#' 
#' For most networks, the default options will be fine, but there is a lot of room
#' for modification if changes are desired
#' 
#' @param data.obj A \code{\link{Cape}} object
#' @param p.or.q The maximum p value (or q value if FDR was used) for significant 
#' main effects and interactions.
#' @param collapsed.net A logical value indicating whether to show the network
#' condensed into linkage blocks (TRUE) or each individual marker (FALSE)
#' @param main A title for the plot
#' @param color.scheme either "DO/CC" or "other". "DO/CC" uses the official "DO/CC"
#' colors for the DO/CC alleles  \url{https://compgen.unc.edu/wp/?page_id=577}
#' "other" uses an unrelated color pallette for multiple alleles.
#' @param pos.col The color to use for positive main effects and interactions
#' must be one of "green", "purple", "red", "orange", "blue", "brown", "yellow", "gray"
#' see \code{\link{get.color}}
#' @param neg.col The color to use for negative main effects and interactions
#' takes the same values as pos.col.
#' @param bg.col The color to be used in pie charts for non-significant main effects.
#' Takes the same values as pos.col
#' @param light.dark Indicates whether pos.col, neg.col, and bg.col should be selected
#' from light colors ("l"), dark colors ("d") or the full spectrum from light to dark ("f")
#' @param node.border.lwd The thickness of the lines around the pie charts
#' @param layout.matrix Users have the option of providing their own layout matrix for the
#' network. This should be a two column matrix indicating the x and y coordinates of each 
#' node in the network.
#' @param zoom Allows the user to zoom in and out on the image if the network is either 
#' running off the edges of the plot or too small in the middle of the plot.
#' @param xshift A constant by which to shift the x values of all nodes in the network.
#' @param yshift A constant by which to shift the y values of all nodes in the network.
#' @param node.radius The size of the pie chart for each node.
#' @param label.nodes A logical value indicating whether the nodes should be labeled.
#' Users may want to remove labels for large networks.
#' @param label.offset The amount by which to offset the node labels from the center of
#' the nodes.
#' @param label.cex The size of the node labels
#' @param legend.radius The size of the legend indicating which pie piece corresponds to which
#' traits.
#' @param legend.cex The size of the labels in the legend.
#' @param legend.position The position of the legend on the plot
#' @param arrow.offset The distance from the center of the node to the arrow head.
#' @param arrow.length The length of the head of the arrow
#' @param edge.lwd The thickness of the arrows showing the interactions.
#'
#' @return This function invisibly returns a list of length two. The first element
#' contains the igraph network object. The second contains the layout matrix for the
#' network. This can be passed in as an argument ("layout.matrix") which provides more
#' control to the user in the layout. Other network layouts from igraph can also be passed 
#' in here.
#'
#' @references Csardi G, Nepusz T: The igraph software package for complex network 
#' research, InterJournal, Complex Systems 1695. 2006. http://igraph.org 
#' 
#' @export
 
plotFullNetwork <- function(data.obj, p.or.q = 0.05,  collapsed.net = TRUE, main = NULL, 
                         color.scheme = c("DO/CC", "other"), pos.col = "brown", neg.col = "blue", 
                         bg.col = "gray", light.dark = "f", node.border.lwd = 1, layout.matrix = NULL, 
                         zoom = 1, xshift = 0, yshift = 0, node.radius = 1, label.nodes = TRUE, 
                         label.offset = 0, label.cex = 1, legend.radius = 1, legend.cex = 1, 
                         legend.position = "topleft", arrow.offset = node.radius, arrow.length = 0.2, 
                         edge.lwd = 2){
  
    
  pheno.tables <- data.obj$max_var_to_pheno_influence
  phenotypes <- names(pheno.tables)	
  
  
  plot.net.edges <- function(net, net.layout, lwd = 1, edge.col = "gray", arrow.offset = 1){
    
    edge.list <- get.edgelist(net, names = FALSE)
    
    if(length(col) == 1){
      edge.col = rep(edge.col, dim(edge.list)[1])
    }
    
    
    for(i in 1:length(edge.list[,1])){
      # cat(i, "\n")
      xy.start = c(net.layout[edge.list[i,1],1], net.layout[edge.list[i,1],2])
      xy.end <- c(net.layout[edge.list[i,2],1], net.layout[edge.list[i,2],2])
      alpha <- 0
      final.xy <- alpha*xy.start + (1-alpha)*(xy.end)
      dist.to.center <- dist(matrix(c(xy.end, final.xy), nrow = 2, byrow = TRUE), method = "euclidean")
      counter <- 0
      while(dist.to.center < arrow.offset && counter < 100){
        alpha <- alpha + 0.01
        final.xy <- alpha*xy.start + (1-alpha)*(xy.end)
        dist.to.center <- dist(matrix(c(xy.end, final.xy), nrow = 2, byrow = TRUE), method = "euclidean")	
        counter = counter + 1
        # print(counter)
      }
      # print(paste("\t", counter))
      # cat(round(xy.start, 2), " : ", round(final.xy, 2), "\n")
      arrows(x0 = xy.start[1], x1 = final.xy[1], y0 = xy.start[2], y1 = final.xy[2], lwd = lwd, col = edge.col[i], length = arrow.length)
    }
    
  }
  
  #This function plots nodes of a graph as polygons where each polygon
  #is subdivided into polygons each with a different color
  #This is for adding information to nodes about which phenotypes they
  #have main effects on
  plot.net.point <- function(x, y, node.radius = node.radius, 
  cols = c("green", "green", "red"), node.label = NULL, label.offset = 0, 
  label.cex = 1, border.col){
    
    draw.pie(x = x, y = y, radius = node.radius, cols = cols, add = TRUE, 
    node.border.lwd = node.border.lwd, border.col = border.col)
    
    
    if(!is.null(node.label)){
      text.x <- x+label.offset; text.y = y+label.offset
      text(x = text.x, y = text.y, labels = node.label, cex = label.cex)
    }
    
  }
  
  
  
  if(collapsed.net){
    blocks <- data.obj$linkage_blocks_collapsed
    net.data <- unlist(data.obj$collapsed_net)
    #convert this into an edge list to be compatible with the uncollapsed network
    sig.locale <- which(net.data != 0, arr.ind = TRUE)
    
    if(length(sig.locale) == 0){
      plot.new()
      plot.window(xlim = c(0,1), ylim = c(0,1))
      text(0.5, 0.5, "No Significant Interactions")
      invisible()
    }
    edge.vals <- net.data[which(net.data != 0)]
    sig.edges <- cbind(sig.locale, edge.vals)
    colnames(sig.edges) <- c("Source", "Target", "Effect")
    block.names <- matrix(NA, ncol = 2, nrow = dim(sig.edges)[1])
    block.names[,1] <- rownames(net.data)[sig.edges[,1]]
    block.names[,2] <- colnames(net.data)[sig.edges[,2]]
    sig.edges[,1:2] <- block.names
    
  }else{
    blocks = data.obj$linkage_blocks_full
    net.data <- data.obj$var_to_var_p_val
    
    var.sig.col <- which(colnames(data.obj$var_to_var_p_val) == "p.adjusted")
    sig.locale <- which(net.data[,var.sig.col] <= p.or.q)
    
    if(length(sig.locale) == 0){
      plot.new()
      plot.window(xlim = c(0,1), ylim = c(0,1))
      text(0.5, 0.5, "No Significant Interactions")
      invisible()
    }
    
    #get the edges for variant to variant interactions
    sig.edges <- matrix(net.data[sig.locale,1:2], nrow = length(sig.locale))
    effect.size <- (as.numeric(net.data[sig.locale,"Effect"])/as.numeric(net.data[sig.locale,"SE"]))
    sig.edges <- cbind(sig.edges, effect.size)
    colnames(sig.edges) <- c("Source", "Target", "Effect")
    
    #add any edges from the phenotype direct influences
    for(ph in 1:length(phenotypes)){
      pheno.stats <- pheno.tables[[ph]][,which(colnames(pheno.tables[[ph]]) != "t.stat")]
      p.adj.locale <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(pheno.stats))))
      sig.ph.edges <- which(as.numeric(pheno.stats[,p.adj.locale]) <= p.or.q)
      if(length(sig.ph.edges) > 0){
        for(r in sig.ph.edges){
          m.coef <- as.numeric(pheno.stats[r,which(colnames(pheno.stats) == "coef")])
          m.se <- as.numeric(pheno.stats[r,which(colnames(pheno.stats) == "se")])
          table.row <- c(pheno.stats[r,1], names(pheno.tables)[ph], m.coef/m.se)
          sig.edges <- rbind(sig.edges, table.row)
        }
      }
    } #end adding phenotype edges
    
  }
  
  
  get.node.cols <- function(node.name, phenotypes, main.effects){
    node.locale <- which(main.effects[,"Source"] == node.name)
    no.effect.col <- bg.col
    node.cols <- rep(no.effect.col, length(phenotypes))
    
    node.effects <- main.effects[node.locale,c("Target", "Effect"),drop=FALSE]
    
    if(length(node.effects) > 0){
      for(i in 1:length(phenotypes)){
        pheno.locale <- which(node.effects[,"Target"] == phenotypes[i])
        if(length(pheno.locale) > 0){
          if(as.numeric(node.effects[pheno.locale,"Effect"]) > 0){node.cols[i] <- get.color(pos.col, light.dark)[2]}
          if(as.numeric(node.effects[pheno.locale,"Effect"]) < 0){node.cols[i] <- get.color(neg.col, light.dark)[2]}
        }
      }
    }
    return(node.cols)
  }
  
  
  main.effects.locale <- which(sig.edges[,2] %in% phenotypes)
  if(length(main.effects.locale) > 0){
    main.effects <- sig.edges[main.effects.locale,,drop=FALSE]
    interactions <- sig.edges[-main.effects.locale,,drop=FALSE]
  }else{
    main.effects <- NULL
    interactions <- sig.edges
  }
  
  
  if(length(interactions) == 0){
    plot.new()
    plot.window(xlim = c(0,1), ylim = c(0,1))
    text(0.5, 0.5, "No Significant Interactions")
    invisible()
  }else{	
    
    edgelist <- matrix(c(as.vector(interactions[,1]), as.vector(interactions[,2])), ncol = 2, byrow = FALSE)
        
    net <- graph.edgelist(edgelist)
    vertex.names <- V(net)$name
    
    if(!label.nodes){vertex.names <- NULL}
    
    edge.col <- rep(NA, ecount(net))
    edge.col[which(as.numeric(interactions[,"Effect"]) < 0)] <- get.color(neg.col, light.dark)[2]
    edge.col[which(as.numeric(interactions[,"Effect"]) > 0)] <- get.color(pos.col, light.dark)[2]
    
    if(is.character(layout.matrix)){
      if(layout.matrix == "manual"){
        tkp.id <- tkplot(net)
        readline(prompt = "Press return when ready:\n")
        layout.matrix <- tkplot.getcoords(tkp.id)
        tkplot.close(tkp.id)
      }else{
        layout.call <- call(layout.matrix, net)
        coord.matrix <- eval(layout.call)
      }			
    }else{ #otherwise layout.matrix is null or a matrix
      if(length(layout.matrix) == 0){
        coord.matrix <- layout.auto(net)
      }else{
        coord.matrix <- layout.matrix	
      }
    }
    #normalize the layout to be between -10 and 10
    norm.coord <- apply(coord.matrix, 2, function(x) x*10/max(abs(x)))
    
    #center the coordinate matrix, so we can spread out the nodes using zoom
    centered.coord <- apply(norm.coord, 2, function(x) x - mean(x))
    # minx <- min(centered.coord[,1]); maxx <- max(centered.coord[,1])
    # miny <- min(centered.coord[,2]); maxy <- max(centered.coord[,2])
    minx <- min(centered.coord); maxx <- max(centered.coord)
    miny <- min(centered.coord); maxy <- max(centered.coord)
    
    
    zoomed.coord <- centered.coord*zoom
    zoomed.coord[,1] <- zoomed.coord[,1] + xshift
    zoomed.coord[,2] <- zoomed.coord[,2] + yshift
    
    block.alleles <- sapply(V(net)$name, function(x) get.block.allele(data.obj, x, collapsed.net))
    allele.cols <- get.allele.colors(color.scheme, unique(block.alleles))
    block.cols <- allele.cols[match(block.alleles, allele.cols[,2]), 3]
    block.cols[which(is.na(block.cols))] <- "black"
    
    
    # par(mar = c(2,2,2,2))
    plot.new()
    plot.window(xlim = c(minx, maxx), ylim = c(miny, maxy))
    par(xpd = TRUE)
    plot.net.edges(net = net, net.layout = zoomed.coord, lwd = edge.lwd, 
    edge.col = edge.col, arrow.offset = arrow.offset)
    
    for(i in 1:length(V(net)$name)){
      plot.net.point(x = zoomed.coord[i,1], y = zoomed.coord[i,2], 
      node.radius = node.radius, 
      cols = get.node.cols(node.name = V(net)$name[i], phenotypes = phenotypes, main.effects = main.effects), 
      node.label = vertex.names[i], label.offset = label.offset, 
      label.cex = label.cex, border.col = block.cols[i])
    }
    
    if(legend.position == "topleft"){
      l.x <- minx; l.y <- maxy
    }
    if(legend.position == "topright"){
      l.x <- maxx; l.y <- maxy
    }
    if(legend.position == "bottomleft"){
      l.x <- minx; l.y <- miny
    }
    if(legend.position == "bottomright"){
      l.x <- maxx; l.y <- miny
    }
    
    
    draw.pie(x = l.x, y = l.y, radius = legend.radius, cols = rep(bg.col, length(phenotypes)), labels = phenotypes, label.cex = legend.cex, node.border.lwd = node.border.lwd, border.col = "black")
    
    if(!is.null(main)){
      text(x = mean(c(maxx, minx)), y = maxy, labels = main)
    }	
    
    par(xpd = FALSE)
    
    results <- list(net, zoomed.coord)
    names(results) <- c("net", "layout.matrix")
    invisible(results)
  }
  
  
}
