#' Calculate m
#' 
#' This function performs the reparameterization
#' step originally described in
#' Carter, G. W., Hays, M., Sherman, A. & Galitski, T. 
#' Use of pleiotropy to model genetic interactions in a 
#' population. PLoS Genet. 8, e1003010 (2012).
#' 
#' @param markers A vector of length two holding the
#' names of the markers being tested.
#' @param beta_m The matrix of beta values from the 
#' pairwise linear regression testing the interaction
#' between the markers.
#'
#' @return This function returns a matrix with one 
#' row holding the m12 and m21 values as described 
#' in Carter et al. 2012
#' 
#' @import corpcor
#' @importFrom stats D
#' @keywords internal

calc_m <- function(markers,beta_m) {
  beta_main <- beta_m[,1:2]
  beta_inter <- beta_m[,3]
  n_rows <- dim(beta_main)[1] #No of ETs
  n_cols <- dim(beta_main)[2]
  
  if(n_rows == n_cols){
    act_delta <- solve(beta_main)%*%beta_inter
  }else{
    tolerance = max(dim(beta_main))*max(D)*.Machine$double.eps
    act_delta <- try(pseudoinverse(beta_main, tol = tolerance)%*%beta_inter, silent = TRUE)
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