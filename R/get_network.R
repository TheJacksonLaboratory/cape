#' Convert the final results to an adjacency 
#' matrix.
#' 
#' This function converts the significant cape 
#' interactions to an adjacency matrix, which 
#' is then used by \code{\link{plot_network}}
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' @param p_or_q A threshold indicating the maximum adjusted p value considered 
#' significant. If an fdr method has been used to correct for multiple testing, 
#' this value specifies the maximum q value considered significant.
#' @param min_std_effect This numerical value offers an additional filtering
#' method. If specified, only interactions with standardized effect sizes greater
#' then the min_std_effect will be returned. 
#' @param standardize A logical value indicating whether the values returned in
#' the adjacency matrix should be effect sizes (FALSE) or standardized effect
#' sizes (TRUE). Defaults to FALSE.
#' @param collapse_linked_markers A logical value. If TRUE markers are combined 
#' into linkage blocks based on correlation. If FALSE, each marker is treated as 
#' an independent observation.
#' @param threshold_power A numerical value indicating the power to which to 
#' raise the marker correlation matrix. This parameter is used in 
#' \code{\link{linkage_blocks_network}} to determine soft thresholding
#' in determining linkage block structure. 
#' Larger values result in more splitting of linkage blocks. Smaller values 
#' result in less splitting. The default value of 1 uses the unmodified
#' correlation matrix to determine linkage block structure.
#' @param verbose A logical value indicating whether to print algorithm progress
#' to standard out.
#' @param plot_linkage_blocks A logical value indicating whether to plot heatmaps
#' showing the marker correlation structure and where the linkage block boundaries
#' were drawn.
#' 
#' @return This function returns the data object with an adjacency matrix defining
#' the final cape network based on the above parameters. The network is put into 
#' the slot collapsed_net if collapse_linked_markers is set to TRUE, and full_net
#' if collapse_linked_markers is set to FALSE. \code{\link{run_cape}} automatically
#' requests both networks be generated.
#' 
#' @examples 
#' \dontrun{
#' data_obj <- get_network(data_obj, geno_obj)
#' }
#' 
#' @import igraph
#' 
#' @export
get_network <- function(data_obj, geno_obj, p_or_q = 0.05, min_std_effect = 0, standardize = FALSE, 
                        collapse_linked_markers = TRUE, threshold_power = 1, verbose = FALSE, 
                        plot_linkage_blocks = FALSE){
  
  if(verbose){cat("Calculating linkage blocks...\n")}
  #get the linkage blocks based on the significant markers
  data_obj <- linkage_blocks_network(data_obj, geno_obj,
    collapse_linked_markers = collapse_linked_markers, threshold_power = threshold_power, 
    plot_blocks = plot_linkage_blocks)
  
  if(collapse_linked_markers){
    blocks <- data_obj$linkage_blocks_collapsed
  }else{
    blocks <- data_obj$linkage_blocks_full
  }
  
  if(length(blocks) == 1){
    stop("There is only one linkage block at this r2 threshold.")
  }
  
  #make sure all covariates are in linkage blocks
  covar_info <- get_covar(data_obj)
  non_allelic_covar <- covar_info$covar_names 
  if(length(non_allelic_covar) > 0){
    for(i in 1:length(non_allelic_covar)){
      covar_locale <- which(names(blocks) == non_allelic_covar[i])
      if(length(covar_locale) == 0){
        blocks[[length(blocks)+1]] <- non_allelic_covar[i]
        names(blocks)[length(blocks)] <- non_allelic_covar[i]
      }	
    }
  }
  
  if(collapse_linked_markers){
    data_obj$linkage_blocks_collapsed <- blocks
  }else{
    data_obj$linkage_blocks_full <- blocks
  }
  
  #build a new network based on the block structure
  all_net_data <- data_obj$var_to_var_p_val
  
  pheno_tables <- data_obj$max_var_to_pheno_influence
  for(i in 1:length(pheno_tables)){
    not_finite <- which(!is.finite(as.numeric(pheno_tables[[i]][,7])))
    pheno_tables[[i]][not_finite,7] <- 0
  }
  phenotypes <- names(pheno_tables)	
  
  rename_with_blocks <- function(marker_names, blocks){
    marker_blocks <- rep(NA, length(marker_names))
    for(i in 1:length(blocks)){
      block_marker_locale <- which(marker_names %in% blocks[[i]])
      marker_blocks[block_marker_locale] <- names(blocks)[i]
    }
    return(marker_blocks)
  }
  

  if(verbose){cat("Creating adjacency matrix...\n")}
  
  #replace each marker name with the name of its block
  source_markers <- rename_with_blocks(marker_names = all_net_data[,1], blocks)
  target_markers <- rename_with_blocks(all_net_data[,2], blocks)
  
  edgelist <- cbind(source_markers, target_markers)
  
  # the igraph::graph_from_edgelist method fails on NA values. remove them.
  not_na <- which(apply(edgelist, 1, function(x) all(!is.na(x))))
  filtered_edgelist <- edgelist[not_na,]
  #filtered_edgelist <- subset(edgelist, !is.na(edgelist[,"source_markers"]) & !is.na(edgelist[,"target_markers"]))
  
  net <- graph_from_edgelist(filtered_edgelist)
  if(standardize){
    weights <- as.numeric(all_net_data[not_na,5])
  }else{
    weights <- as.numeric(all_net_data[not_na,3])
  }
  var_sig_col <- 7
  non_sig_locale <- which(as.numeric(all_net_data[not_na,var_sig_col]) > p_or_q)
  weights[non_sig_locale] <- 0
  E(net)$weight <- weights
  
  #remove multiple edges between blocks, take the edge with the maximum magnitude
  edge_attr_comb = list(weight = function(x) x[which.max(abs(x))], name="ignore")
  simple_net <- simplify(net, remove.multiple = TRUE, remove.loops = FALSE, edge.attr.comb = edge_attr_comb)
  
  adj_mat <- as.matrix(as_adjacency_matrix(simple_net, type = "both", attr = "weight"))
  
  #put the adjacency matrix in chromosome order
  split_labels <- unlist(lapply(strsplit(rownames(adj_mat), "Chr"), function(x) x[2]))
  split_again <- strsplit(split_labels, "_")
  chr_num <- unlist(lapply(split_again, function(x) x[1]))
  block_num <- unlist(lapply(split_again, function(x) x[2]))
  chr_order <- match(chr_num, unique(data_obj$chromosome))
  new_order <- sort_by_then_by(tableX = cbind(chr_order, block_num), sort_cols = c(1,2), col_type = c("n", "n"), return_order = TRUE)
  
  #put in chromosome order
  adj_mat <- adj_mat[new_order[,1],]
  adj_mat <- adj_mat[new_order[,2],]
  adj_mat <- adj_mat[,new_order[,1]]
  adj_mat <- adj_mat[,new_order[,2]]
  
  
  if(verbose){cat("Adding main effects...\n")}
  
  #Now add the phenotypic effects continuing to use the maximum significant effect from each block
  pheno_mat <- matrix(0, nrow = length(blocks), ncol = length(phenotypes))
  colnames(pheno_mat) <- phenotypes
  rownames(pheno_mat) <- names(blocks)
  
  pheno_sig_col <- which(colnames(pheno_tables[[1]]) == "p_adjusted")
  
  get_block_inf <- function(block){
    all_markers <- blocks[[block]]
    
    for(i in 1:length(pheno_tables)){
      sig_inf <- pheno_tables[[i]][which(as.numeric(pheno_tables[[i]][,pheno_sig_col]) <= p_or_q),,drop = FALSE]
      if(length(sig_inf) > 0){
        block_locale <- which(sig_inf[,1] %in% all_markers)
        if(length(block_locale) > 0){
          all_effects <- as.numeric(sig_inf[block_locale,"coef"])/as.numeric(sig_inf[block_locale,"se"])
          pheno_mat[block,i] <- all_effects[which(abs(all_effects) == max(abs(all_effects)))][1]
        }
      }
    }
    return(pheno_mat)	
  }
  
  for(i in 1:length(blocks)){
    pheno_mat <- get_block_inf(block = names(blocks)[i])
  }
  
  #some alleles are never tested because they do not have any variance
  #but they still get a slot in the pheno_mat. Remove these before
  #binding the results together
  were_tested <- intersect(rownames(pheno_mat), rownames(adj_mat))
  tested_locale <- match(were_tested, rownames(pheno_mat))
  pheno_mat <- pheno_mat[tested_locale,]
  
  final_mat <- cbind(adj_mat, pheno_mat)
  
  if(collapse_linked_markers){
    data_obj$collapsed_net <- final_mat
  }else{
    data_obj$full_net <- final_mat	
  }
  
  return(data_obj)
  
}
