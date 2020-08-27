#' Check selected markers for linear independence.
#' 
#' This function checks a marker matrix selected by
#' \code{\link{select.markers.for.pairscan}} for linear
#' independence by checking the correlation between
#' pairs of markers. If any are perfectly correlated,
#' only one marker from the block of perfect correlation
#' is kept.
#'
#' @param data.obj a \code{\link{Cape}} object
#'
#' @return This function returns a list with two elements.
#' independent.markers: A matrix of the markers that are 
#' linearly independent. 
#' rejected.markers: A vector indicating which markers were
#' removed for violating linear independence.
#' 

get.linearly.independent <- function(data.obj){
  
  matrixX <- data.obj$geno_for_pairscan
  
  if(dim(matrixX)[2] == 1){
    return(matrixX)
  }
  
  #use precision to the 3rd decimal place
  matrixX <- round(matrixX, 3)
  
  #find the markers without variation
  num.geno <- apply(matrixX, 2, function(x) length(unique(x[!is.na(x)])))
  bad.markers <- which(num.geno < 2)
  rejected.markers <- names(num.geno[num.geno < 2])
  
  if(length(rejected.markers) > 0){
    matrixX <- matrixX[,-bad.markers]
    rejected.markers <- bad.markers
  }
  
  all.cor <- cor(matrixX)
  all.cor[lower.tri(all.cor, diag = TRUE)] <- 0
  # hist(all.cor)
  
  
  #find blocks of perfect correlation, and keep one
  #marker from each block
  
  perfect.cor <- apply(all.cor, 1, function(x) which(x == 1))
  
  #define blocks 
  if(length(perfect.cor) > 0){
    cor.blocks <- list()
    if(length(perfect.cor[[1]]) > 0){
      cor.blocks[[1]] <- names(perfect.cor[[1]])
    }else{
      cor.blocks[[1]] <- names(perfect.cor)[1]	
    }
    block.num <- 1
    for(i in 2:length(perfect.cor)){
      common.markers <- length(intersect(names(perfect.cor[[i]]), names(perfect.cor[[(i-1)]])))
      if(common.markers > 0){
        cor.blocks[[block.num]] <- unique(c(cor.blocks[[block.num]], names(perfect.cor[[i]])))
      }else{
        block.num <- block.num + 1
        cor.blocks[[block.num]] <- names(perfect.cor)[i]	
      }
    }
    
    #go through the blocks and take the first element of each
    uncor.markers <- unlist(lapply(cor.blocks, function(x) x[1]))
    new.matrix <- matrixX[,uncor.markers]
    
    rejected.markers <- setdiff(colnames(matrixX), uncor.markers)
  }else{
    new.matrix <- matrixX
    rejected.markers <- NULL	
  }
  
  results <- list(new.matrix, rejected.markers)
  names(results) <- c("independent.markers", "rejected.markers")
  return(results)	
  
}