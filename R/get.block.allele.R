#' Get allele assignments for linkage blocks
#'
#' This is an internal function used by \link{\code{plotFullNetwork}}
#' to retrieve the allele names contained in linkage blocks for 
#' plotting using the appropriate colors.
#' 
#' @param data.obj The cape object
#' @param block.name The name of the linkage block as defined
#' by the name of the list called linkage_blocks_full or 
#' linkage_blocks_collapsed
#' @param collapsed.net A logical value indicating whether the 
#' block names to be used are from linkage_blocks_collapsed (TRUE)
#' or linkage_blocks_full (FALSE)
#' 
#' @return A vector of the alleles contained in the linkage block

get.block.allele <- function(data.obj, block.name, collapsed.net){
  
  if(collapsed.net){
    blocks <- data.obj$linkage_blocks_collapsed
  }else{
    blocks <- data.obj$linkage_blocks_full	
  }
  block.locale <- which(names(blocks) == block.name)
  block.allele <- strsplit(blocks[[block.locale]], "_")[[1]][2]
  
  return(block.allele)
}
