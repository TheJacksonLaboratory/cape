#' Convert the final results to a form plotted by \code{\link{plotNetwork}} and \code{\link{plotCollapsedVarInf}}
#' 
#' This function collapses a network based on linkage.
#' It returns the data object with a new weighted 
#' adjacency matrix in which the weights are 
#' standardized effects of each influence. If there
#' is more than one marker in a block, the weight
#' is the maximum influence between markers in each
#' block
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param p.or.q A threshold indicating the maximum adjusted p value considered 
#' significant. If an fdr method has been used to correct for multiple testing, 
#' this value specifies the maximum q value considered significant.
#' @param min.std.effect 
#' @param standardize A logical value. If FALSE, the interaction terms are stored 
#' as marker to marker influence coefficients. If TRUE, the coefficients are 
#' standardized by their standard errors.
#' @param collapse.linked.markers A logical value. If TRUE markers are combined 
#' into linkage blocks based on correlation. If FALSE, each marker is treated as 
#' an independent observation.
#' @param threshold.power 
#' @param verbose default = FALSE
#' @param plot.linkage.blocks 
#' @param lookup.marker.position 
#' 
get.network <- function(data.obj, p.or.q = 0.05, min.std.effect = 0, standardize = FALSE, 
                        collapse.linked.markers = TRUE, threshold.power = 1, verbose = FALSE, 
                        plot.linkage.blocks = FALSE, lookup.marker.position = FALSE){
  
  if(verbose){cat("Calculating linkage blocks...\n")}
  #get the linkage blocks based on the significant markers
  data.obj <- linkage.blocks.network(data.obj, collapse.linked.markers = collapse.linked.markers, 
                                     threshold.power = threshold.power, plot.blocks = plot.linkage.blocks, 
                                     lookup.marker.position = lookup.marker.position)
  
  if(collapse.linked.markers){
    blocks <- data.obj$linkage_blocks_collapsed
  }else{
    blocks <- data.obj$linkage_blocks_full
  }
  
  if(length(blocks) == 1){
    stop("There is only one linkage block at this r2 threshold.")
  }
  
  #make sure all covariates are in linkage blocks
  covar.info <- get.covar(data.obj)
  non.allelic.covar <- covar.info$covar.names 
  if(length(non.allelic.covar) > 0){
    for(i in 1:length(non.allelic.covar)){
      covar.locale <- which(names(blocks) == non.allelic.covar[i])
      if(length(covar.locale) == 0){
        blocks[[length(blocks)+1]] <- non.allelic.covar[i]
        names(blocks)[length(blocks)] <- non.allelic.covar[i]
      }	
    }
  }
  
  if(collapse.linked.markers){
    data.obj$linkage_blocks_collapsed <- blocks
  }else{
    data.obj$linkage_blocks_full <- blocks
  }
  
  #build a new network based on the block structure
  all.net.data <- data.obj$var_to_var_p_val
  
  pheno.tables <- data.obj$max_var_to_pheno_influence
  for(i in 1:length(pheno.tables)){
    not.finite <- which(!is.finite(as.numeric(pheno.tables[[i]][,7])))
    pheno.tables[[i]][not.finite,7] <- 0
  }
  phenotypes <- names(pheno.tables)	
  
  rename.with.blocks <- function(marker.names, blocks){
    marker.blocks <- rep(NA, length(marker.names))
    for(i in 1:length(blocks)){
      block.marker.locale <- which(marker.names %in% blocks[[i]])
      marker.blocks[block.marker.locale] <- names(blocks)[i]
    }
    return(marker.blocks)
  }
  
  if(verbose){cat("Creating adjacency matrix...\n")}
  
  #replace each marker name with the name of its block
  source.markers <- rename.with.blocks(all.net.data[,1], blocks)
  target.markers <- rename.with.blocks(all.net.data[,2], blocks)
  
  edgelist <- cbind(source.markers, target.markers)
  
  # the igraph::graph_from_edgelist method fails on NA values. remove them.
  filtered_edgelist <- subset(edgelist, !is.na(edgelist[,"source.markers"]) & !is.na(edgelist[,"target.markers"]))
  
  net <- igraph::graph_from_edgelist(filtered_edgelist)
  if(standardize){
    weights <- as.numeric(all.net.data[,5])
  }else{
    weights <- as.numeric(all.net.data[,3])
  }
  var.sig.col <- 7
  non.sig.locale <- which(as.numeric(all.net.data[,var.sig.col]) > p.or.q)
  weights[non.sig.locale] <- 0
  igraph::E(net)$weight <- weights
  
  #remove multiple edges between blocks, take the edge with the maximum magnitude
  edge.attr.comb = list(weight = function(x) x[which.max(abs(x))], name="ignore")
  simple.net <- igraph::simplify(net, remove.multiple = TRUE, remove.loops = FALSE, edge.attr.comb = edge.attr.comb)
  
  adj.mat <- as.matrix(igraph::as_adjacency_matrix(simple.net, type = "both", attr = "weight"))
  
  #put the adjacency matrix in chromosome order
  split.labels <- unlist(lapply(strsplit(rownames(adj.mat), "Chr"), function(x) x[2]))
  split.again <- strsplit(split.labels, "_")
  chr.num <- unlist(lapply(split.again, function(x) x[1]))
  block.num <- unlist(lapply(split.again, function(x) x[2]))
  chr.order <- match(chr.num, unique(data.obj$chromosome))
  new.order <- sortByThenBy(tableX = cbind(chr.order, block.num), sort.cols = c(1,2), col.type = c("n", "n"), return.order = TRUE)
  
  #put in chromosome order
  adj.mat <- adj.mat[new.order[,1],]
  adj.mat <- adj.mat[new.order[,2],]
  adj.mat <- adj.mat[,new.order[,1]]
  adj.mat <- adj.mat[,new.order[,2]]
  
  
  if(verbose){cat("Adding main effects...\n")}
  
  #Now add the phenotypic effects continuing to use the maximum significant effect from each block
  pheno.mat <- matrix(0, nrow = length(blocks), ncol = length(phenotypes))
  colnames(pheno.mat) <- phenotypes
  rownames(pheno.mat) <- names(blocks)
  
  pheno.sig.col <- which(colnames(pheno.tables[[1]]) == "p.adjusted")
  
  get.block.inf <- function(block){
    all.markers <- blocks[[block]]
    
    for(i in 1:length(pheno.tables)){
      sig.inf <- pheno.tables[[i]][which(as.numeric(pheno.tables[[i]][,pheno.sig.col]) <= p.or.q),,drop = FALSE]
      if(length(sig.inf) > 0){
        block.locale <- which(sig.inf[,1] %in% all.markers)
        if(length(block.locale) > 0){
          all.effects <- as.numeric(sig.inf[block.locale,"coef"])/as.numeric(sig.inf[block.locale,"se"])
          pheno.mat[block,i] <- all.effects[which(abs(all.effects) == max(abs(all.effects)))][1]
        }
      }
    }
    return(pheno.mat)	
  }
  
  for(i in 1:length(blocks)){
    pheno.mat <- get.block.inf(block = names(blocks)[i])
  }
  
  #some alleles are never tested because they do not have any variance
  #but they still get a slot in the pheno.mat. Remove these before
  #binding the results together
  not.tested <- setdiff(rownames(pheno.mat), rownames(adj.mat))
  if(length(not.tested) > 0){
    not.tested.locale <- match(not.tested, rownames(pheno.mat))
    pheno.mat <- pheno.mat[-not.tested.locale,]
  }
  
  final.mat <- cbind(adj.mat, pheno.mat)
  
  if(collapse.linked.markers){
    data.obj$collapsed_net <- final.mat
  }else{
    data.obj$full_net <- final.mat	
  }
  
  return(data.obj)
  
}






