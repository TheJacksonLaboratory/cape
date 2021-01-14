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
#' @param data_obj a \code{\link{Cape}} object
#' @param traits_which A vector of integers, of at least length two
#'  specifying which eigentraits should be analyzed.
#'
#' @seealso \code{\link{plot_svd}}
#'
#' @return updated \code{\link{Cape}} object
#' 
#' @examples 
#' \dontrun{
#' data_obj <- selecct_eigentraits(data_obj, traits_which = 1:3)
#' }
#'
#' @export
select_eigentraits <- function(data_obj, traits_which = c(1,2)){
  
  check_underscore(data_obj)
  # check_bad_markers(data_obj)
  
  ET <- data_obj$ET
  
  if(length(traits_which) > ncol(ET)){stop("There are more eigentraits specified than exist. Please select a smaller number in the parameter file.")}
  
  selected_ET <- ET[,traits_which]
  data_obj$ET <- selected_ET
  return(data_obj)
  
}