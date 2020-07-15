#' Calculate m
#' 
#' This function performs the reparametrization
#' step originally described in
#' Carter, G. W., Hays, M., Sherman, A. & Galitski, T. 
#' Use of pleiotropy to model genetic interactions in a 
#' population. PLoS Genet. 8, e1003010 (2012).
#' 
#' @param markers A vector of length two holding the
#' names of the markers being tested.
#' @param beta.m The matrix of beta values from the 
#' pairwise linear regression testing the interaction
#' between the markers.
#'
#' @return This function returns a matrix with one 
#' row holding the m12 and m21 values as described 
#' in Carter et al. 2012


calc.m <- function(markers,beta.m) {
  beta.main <- beta.m[,1:2]
  beta.inter <- beta.m[,3]
  n.rows <- dim(beta.main)[1] #No of ETs
  n.cols <- dim(beta.main)[2]
  
  if(n.rows == n.cols){
    act_delta <- solve(beta.main)%*%beta.inter
  }else{
    tolerance = max(dim(beta.main))*max(D)*.Machine$double.eps
    act_delta <- try(corpcor::pseudoinverse(beta.main, tol = tolerance)%*%beta.inter, silent = TRUE)
    if(class(act_delta) == "try-error"){
      act_delta <- c(NA, NA)
    }
  }
  
  m12 = act_delta[1]/(1+act_delta[2])
  m21 = act_delta[2]/(1+act_delta[1])
  
  results <- cbind(markers[1],markers[2],m12,m21)
  colnames(results) <- c("marker 1","marker 2","m12","m21")
  rownames(results)  <- NULL
  return(results)
  
}