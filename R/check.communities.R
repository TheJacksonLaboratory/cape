#' Check community assignments
#' 
#' This function checks community assignments
#' and makes sure that communities are not 
#' split up along the chromosome
#' 
#' @param comm community from \code{fastgreedy.community}
#' 
check.communities <- function(comm){
  
  new.comm <- comm
  u_comm <- unique(comm)
  com.pos <- lapply(u_comm, function(x) which(comm == x))
  com.consec <- lapply(com.pos, consec.pairs)
  com.dist <- lapply(com.consec, function(x) x[,2] - x[,1])
  big.splits <- lapply(com.dist, function(x) which(x > 1))
  
  next.comm <- max(comm) + 1
  for(i in 1:length(big.splits)){
    if(length(big.splits[[i]]) > 0){
      num.splits <- length(big.splits[[i]])
      for(j in 1:length(big.splits[[i]])){
        break.start <- com.pos[[i]][(big.splits[[i]][j]+1)]
        if(j == num.splits){
          break.end <- com.pos[[i]][length(com.pos[[i]])]
        }else{
          break.end <- com.pos[[i]][(big.splits[[i]][(j+1)])]
        }
        new.comm[break.start:break.end] <- next.comm
        #cbind(comm, new.comm)
        next.comm = next.comm + 1
      }
    }
  }
  return(new.comm)	
}