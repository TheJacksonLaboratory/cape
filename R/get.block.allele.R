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
