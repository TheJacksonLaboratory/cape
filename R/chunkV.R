#' Bins a vector into chunks
#' 
#' This function is used to chunk a vector into 
#' bins. It is primarily used to divide jobs 
#' sent to a certain number of parallel clusters. 

#' @param V The vector to be divided into chunks
#' @param num_chunks The number of chunks in which
#' to bin the vector.
#'
#' @return This function returns a list of length 
#' num_chunks. Each element in the list contains
#' the elements of the vector that were assigned 
#' to that bin.
#' @keywords internal

chunkV <- function(V, num_chunks){
  
  v_list <- vector(mode = "list", length = num_chunks)
  num_per_chunk <- floor(length(V)/num_chunks)
  num_leftover <- length(V)%%num_chunks
  
  #distribute the leftovers througout the list
  without_leftovers <- num_chunks - num_leftover
  extra_per_chunk <- c(rep(0, without_leftovers), rep(1, num_leftover))
  
  start_chunk <- 1
  for(i in 1:num_chunks){
    total <- num_per_chunk + extra_per_chunk[i]
    if(i < num_chunks){
      v_list[[i]] <- V[start_chunk:(start_chunk+total-1)]
      start_chunk <- start_chunk+total
    }else{
      v_list[[i]] <- V[start_chunk:length(V)]
    }
  }
  
  
  return(v_list)
  
}