calc.m <- function(markers,beta.m,beta.cov) {
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