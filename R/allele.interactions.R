#' Makes a table of interacting alleles
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param collapsed.net boolean, default: TRUE
#' @param color.scheme string "CC/DO" or "other"
#'
#' @return updated \code{\link{Cape}} object
#'
#' @export
allele.interactions <- function(data.obj, collapsed.net = TRUE, color.scheme = c("CC/DO", "other")){
  
  
  require(igraph)	
  
  if(collapsed.net){
    net <- data.obj$collapsed.net
    pdf.text <- "Collapsed.Net"
  }else{
    net <- data.obj$full.net	
    pdf.text = "Full.Net"
  }
  
  covar.info <- get.covar(data.obj)
  covar.locale <- which(rownames(net) %in% covar.info$covar.names)
  
  just.markers <- net[,1:nrow(net)]
  just.main <- net[,(nrow(net)+1):ncol(net)]
  
  markers <- rownames(just.markers)
  
  
  blocks <- rownames(net)
  alleles <- sapply(markers, function(x) get.block.allele(data.obj, x, collapsed.net))
  na.locale <- which(is.na(alleles))
  alleles[na.locale] <- names(alleles)[na.locale]
  
  u_alleles <- unique(alleles[-na.locale])
  
  allele.colors <- get.allele.colors(color.scheme, u_alleles)
  allele.abbr <- allele.colors[,2]
  
  #============================================================
  #internal functions
  #============================================================
  plot.counts <- function(counts, label){
    
    allele.locale <- match(names(counts), allele.abbr)
    bar.colors <- allele.colors[allele.locale,3]
    bar.colors[which(is.na(bar.colors))] <- "gray"
    bar.order <- order(counts)
    max.count <- max(counts)
    
    a <- barplot(counts[bar.order], col = bar.colors[bar.order], names = names(counts)[bar.order], ylim = c(0, max.count*1.25), cex.names = 1.5, main = label)	
    text(x = a, y = counts[bar.order]+(max.count*0.05), labels = counts[bar.order], cex = 1.5)
  }
  
  plot.table <- function(net.table){
    inter.net <- graph.adjacency(net.table, "directed", weighted = TRUE)
    node.names <- V(inter.net)$name
    allele.locale <- match(node.names, allele.abbr)
    bar.colors <- allele.colors[allele.locale,3]
    bar.colors[which(is.na(bar.colors))] <- "gray"
    V(inter.net)$color <- bar.colors
    plot(inter.net, edge.width = E(inter.net)$weight, layout = layout.star(inter.net))
    invisible(inter.net)
  }
  #============================================================
  
  rownames(just.markers) <- colnames(just.markers) <- alleles
  alleles.as.sources <- apply(just.markers, 1, function(x) length(which(x != 0)))
  alleles.as.sources.counts <- sapply(u_alleles, function(x) sum(alleles.as.sources[which(names(alleles.as.sources) == x)]))
  
  
  alleles.as.targets <- apply(just.markers, 2, function(x) length(which(x != 0)))
  alleles.as.targets.counts <- sapply(u_alleles, function(x) sum(alleles.as.targets[which(names(alleles.as.targets) == x)]))
  
  
  pdf(paste("Allele.Interactions.", pdf.text, ".pdf", sep = ""), width = 8, height = 6)
  # quartz(width = 10, height = 6)
  plot.counts(alleles.as.sources.counts, "Alleles As Sources")
  # quartz(width = 10, height = 6)
  plot.counts(alleles.as.targets.counts, "Alleles As Targets")
  
  #count how many times each allele interacts with an allele from the same strain
  #get a more detailed table of counts for individual interactions
  intra.table <- rep(0, length(u_alleles))
  names(intra.table) <- u_alleles
  for(i in 1:length(u_alleles)){
    allele.locale <- which(alleles == u_alleles[i])
    if(length(allele.locale) > 0){
      allele.allele.interactions <- just.markers[allele.locale,allele.locale,drop=FALSE]
      intra.table[i] <- length(which(allele.allele.interactions != 0))
    }
  }
  plot.counts(intra.table, "Intra-Strain Interactions")
  
  #get a more detailed table of counts for individual interactions
  inter.table <- matrix(0, nrow = (length(u_alleles)+length(covar.names)), ncol = (length(u_alleles)+length(covar.names)))
  rownames(inter.table) <- colnames(inter.table) <- c(u_alleles, covar.names)
  for(i in 1:length(u_alleles)){
    allele.locale <- which(alleles == u_alleles[i])
    if(length(allele.locale) > 0){
      allele.as.source <- just.markers[allele.locale,,drop=FALSE]
      allele.as.target <- just.markers[,allele.locale,drop=FALSE]			
      sig.as.source.locale <- which(allele.as.source != 0, arr.ind = TRUE)
      sig.as.target.locale <- which(allele.as.target != 0, arr.ind = TRUE)
      if(length(sig.as.source.locale) > 0){
        sig.targets <- alleles[sig.as.source.locale[,2]]
        sig.sources <- alleles[sig.as.target.locale[,1]]
        target.counts <- table(sig.targets)
        source.counts <- table(sig.sources)
        inter.table[i,names(target.counts)] <- target.counts
        inter.table[names(source.counts), i] <- source.counts
      }
    }
  }
  
  
  
  net <- plot.table(inter.table)
  undir.net <- as.undirected(net)
  layout.mat <- get.layout.mat(length(u_alleles), "landscape")
  layout(layout.mat)
  par(mar = c(1,1,1,1))
  for(i in 1:length(u_alleles)){
    others <- neighbors(undir.net, u_alleles[i])
    sub.net <- induced_subgraph(net, c(u_alleles[i], names(others)))
    plot(sub.net, main = u_alleles[i], layout = layout.star(sub.net))
  }
  
  
  deg.in <- sort(degree(net, mode = "in"))
  deg.out <- sort(degree(net, mode = "out"))
  deg.all <- sort(degree(net, mode = "all"))
  
  a <- heatmap(inter.table)
  # quartz(width = 8, height = 6)
  layout(matrix(c(1:3), nrow = 1), widths = c(0.2, 1, 0.2))
  par(mar = c(5,4,4,2))
  plot.new();plot.window(xlim = c(0, 1), ylim = c(0, 1))
  imageWithText(inter.table[a$rowInd,a$colInd], col.names = colnames(inter.table)[a$colInd], row.names = rownames(inter.table)[a$rowInd], col.text.rotation = 0, cex = 1, row.text.cex = 1.5, col.text.cex = 1.5, col.text.adj = 0.5, xlab = "Target", ylab = "Source")
  plot.new();plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  clust.table <- hclust(dist(inter.table))
  table.order <- clust.table$order
  layout(matrix(c(1:3), nrow = 1), widths = c(0.2, 1, 0.2))
  par(mar = c(5,4,4,2))
  plot.new();plot.window(xlim = c(0, 1), ylim = c(0, 1))
  imageWithText(inter.table[table.order,table.order], col.names = colnames(inter.table)[table.order], row.names = rownames(inter.table)[table.order], col.text.rotation = 0, cex = 1, row.text.cex = 1.5, col.text.cex = 1.5, col.text.adj = 0.5, xlab = "Target", ylab = "Source", col.scale = "blue")
  plot.new();plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  dev.off()
  
}