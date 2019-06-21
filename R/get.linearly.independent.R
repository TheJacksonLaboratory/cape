#' Find linearly independent markers based on effect size
#' 
#' This function is a specialized script for getting
#' linearly independent markers based on effect size
#' instead of going through a matrix in order, this
#' script builds a matrix of linearly independent markers
#' starting with markers with the highest effect size
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param verbose boolean, default = FALSE
#'
#' @export
get.linearly.independent <- function(data.obj, verbose = FALSE){
  
  matrixX <- data.obj$geno_for_pairscan
  
  if(dim(matrixX)[2] == 1){
    return(matrixX)
  }
  
  #use precision to the 3rd decimal place
  orig.matrixX <- matrixX
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
      if(verbose){report.progress(i, length(perfect.cor))}
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