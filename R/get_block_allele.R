#' Get allele assignments for linkage blocks
#'
#' This is an internal function used by \code{\link{plot_full_network}}
#' to retrieve the allele names contained in linkage blocks for 
#' plotting using the appropriate colors.
#' 
#' @param data_obj The cape object
#' @param block_name The name of the linkage block as defined
#' by the name of the list called linkage_blocks_full or 
#' linkage_blocks_collapsed
#' @param collapsed_net A logical value indicating whether the 
#' block names to be used are from linkage_blocks_collapsed (TRUE)
#' or linkage_blocks_full (FALSE)
#' 
#' @return A vector of the alleles contained in the linkage block
#' @keywords internal

get_block_allele <- function(data_obj, block_name, collapsed_net){
  
  if(collapsed_net){
    blocks <- data_obj$linkage_blocks_collapsed
  }else{
    blocks <- data_obj$linkage_blocks_full	
  }
  block_locale <- which(names(blocks) == block_name)
  block_allele <- strsplit(blocks[[block_locale]], "_")[[1]][2]
  
  return(block_allele)
}
