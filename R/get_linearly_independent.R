#' Check selected markers for linear independence.
#' 
#' This function checks a marker matrix selected by
#' \code{\link{select_markers_for_pairscan}} for linear
#' independence by checking the correlation between
#' pairs of markers. If any are perfectly correlated,
#' only one marker from the block of perfect correlation
#' is kept.
#'
#' @param data_obj a \code{\link{Cape}} object
#'
#' @return This function returns a list with two elements.
#' independent_markers: A matrix of the markers that are 
#' linearly independent. 
#' rejected_markers: A vector indicating which markers were
#' removed for violating linear independence.
#' @keywords internal
#' 

get_linearly_independent <- function(data_obj){
  
  matrixX <- data_obj$geno_for_pairscan
  
  if(dim(matrixX)[2] == 1){
    return(matrixX)
  }
  
  #use precision to the 3rd decimal place
  matrixX <- round(matrixX, 3)
  
  #find the markers without variation
  num_geno <- apply(matrixX, 2, function(x) length(unique(x[!is.na(x)])))
  bad_markers <- which(num_geno < 2)
  rejected_markers <- names(num_geno[num_geno < 2])
  
  if(length(rejected_markers) > 0){
    matrixX <- matrixX[,-bad_markers]
    rejected_markers <- bad_markers
  }
  
  all_cor <- cor(matrixX)
  all_cor[lower.tri(all_cor, diag = TRUE)] <- 0
  # hist(all_cor)
  
  
  #find blocks of perfect correlation, and keep one
  #marker from each block
  
  perfect_cor <- apply(all_cor, 1, function(x) which(x == 1))
  
  #define blocks 
  if(length(perfect_cor) > 0){
    cor_blocks <- list()
    if(length(perfect_cor[[1]]) > 0){
      cor_blocks[[1]] <- names(perfect_cor[[1]])
    }else{
      cor_blocks[[1]] <- names(perfect_cor)[1]	
    }
    block_num <- 1
    for(i in 2:length(perfect_cor)){
      common_markers <- length(intersect(names(perfect_cor[[i]]), names(perfect_cor[[(i-1)]])))
      if(common_markers > 0){
        cor_blocks[[block_num]] <- unique(c(cor_blocks[[block_num]], names(perfect_cor[[i]])))
      }else{
        block_num <- block_num + 1
        cor_blocks[[block_num]] <- names(perfect_cor)[i]	
      }
    }
    
    #go through the blocks and take the first element of each
    uncor_markers <- unlist(lapply(cor_blocks, function(x) x[1]))
    new_matrix <- matrixX[,uncor_markers]
    
    rejected_markers <- setdiff(colnames(matrixX), uncor_markers)
  }else{
    new_matrix <- matrixX
    rejected_markers <- NULL	
  }
  
  results <- list(new_matrix, rejected_markers)
  names(results) <- c("independent_markers", "rejected_markers")
  return(results)	
  
}