#' Assign selected eigentraits in the Cape object
#' 
#' This function is used to identify which eigentraits
#' will be analyzed in the Cape run. After eigentrait 
#' decomposition of n traits, there will be n eigentraits.
#' If there are more than two eigentraits, the user may
#' wish to analyze a subset of them. This function specifies
#' which of the eigentraits will be analyzed by Cape. It does
#' this by subsetting the ET matrix to only those eigentraits
#' specified. The traits not selected are deleted from the object.
#' 
#' @param data.obj a \code{\link{Cape}} object
#' @param traits.which A vector of integers, of at least length two
#'  specifying which eigentraits should be analyzed.
#'
#' @seealso \link{\code{plotSVD}}
#'
#' @return updated \code{\link{Cape}} object
#'
#' @export
#' 
select.eigentraits <- function(data.obj, traits.which = c(1,2)){
  
  check.underscore(data.obj)
  # check.bad.markers(data.obj)
  
  ET <- data.obj$ET
  
  if(length(traits.which) > ncol(ET)){stop("There are more eigentraits specified than exist. Please select a smaller number in the parameter file.")}
  
  selected.ET <- ET[,traits.which]
  data.obj$ET <- selected.ET
  return(data.obj)
  
}