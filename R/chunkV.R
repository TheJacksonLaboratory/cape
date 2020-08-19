#' Bins a vector into chunks
#' 
#' This function is used to chunk a vector into 
#' bins. It is primarily used to divide jobs 
#' sent to a certain number of parallel clusters. 

#' @param V The vector to be divided into chunks
#' @param num.chunks The number of chunks in which
#' to bin the vector.
#'
#' @return This function returns a list of length 
#' num.chunks. Each element in the list contains
#' the elements of the vector that were assigned 
#' to that bin.

chunkV <- function(V, num.chunks){
  
  v.list <- vector(mode = "list", length = num.chunks)
  num.per.chunk <- floor(length(V)/num.chunks)
  num.leftover <- length(V)%%num.chunks
  
  #distribute the leftovers througout the list
  without.leftovers <- num.chunks - num.leftover
  extra.per.chunk <- c(rep(0, without.leftovers), rep(1, num.leftover))
  
  start.chunk <- 1
  for(i in 1:num.chunks){
    total <- num.per.chunk + extra.per.chunk[i]
    if(i < num.chunks){
      v.list[[i]] <- V[start.chunk:(start.chunk+total-1)]
      start.chunk <- start.chunk+total
    }else{
      v.list[[i]] <- V[start.chunk:length(V)]
    }
  }
  
  
  return(v.list)
  
}