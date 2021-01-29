#' Check community assignments
#' 
#' This function checks community assignments
#' and makes sure that communities are not 
#' split up along the chromosome
#' 
#' @param comm community from \code{fastgreedy.community}
#' @keywords internal
#' 
check_communities <- function(comm){
  
  new_comm <- comm
  u_comm <- unique(comm)
  com_pos <- lapply(u_comm, function(x) which(comm == x))
  com_consec <- lapply(com_pos, consec_pairs)
  com_dist <- lapply(com_consec, function(x) x[,2] - x[,1])
  big_splits <- lapply(com_dist, function(x) which(x > 1))
  
  next_comm <- max(comm) + 1
  for(i in 1:length(big_splits)){
    if(length(big_splits[[i]]) > 0){
      num_splits <- length(big_splits[[i]])
      for(j in 1:length(big_splits[[i]])){
        break_start <- com_pos[[i]][(big_splits[[i]][j]+1)]
        if(j == num_splits){
          break_end <- com_pos[[i]][length(com_pos[[i]])]
        }else{
          break_end <- com_pos[[i]][(big_splits[[i]][(j+1)])]
        }
        new_comm[break_start:break_end] <- next_comm
        #cbind(comm, new_comm)
        next_comm = next_comm + 1
      }
    }
  }
  return(new_comm)	
}