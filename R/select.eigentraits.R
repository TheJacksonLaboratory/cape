#' Assign selected eigentraits in the Cape object
#' 
#' This script is used to select individual eigentraits
#' after viewing the results of the svd
#' The script defaults to eigentraits 1 and 2
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param traits.which integer or array of integers selecting the eigentrait column
#'
#' @return updated \code{\link{Cape}} object
select.eigentraits <- function(data.obj, traits.which = c(1,2)){
  
  check.underscore(data.obj)
  # check.bad.markers(data.obj)
  
  ET <- data.obj$ET
  
  if(length(traits.which) > ncol(ET)){stop("There are more eigentraits specified than exist. Please select a smaller number in the parameter file.")}
  
  selected.ET <- ET[,traits.which]
  data.obj$ET <- selected.ET
  return(data.obj)
  
}