#' This function bins a vector into a given number of chunks
#' 
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