#' pheatmap generate breaks found at this link
#' https://cran.r-project.org/package=pheatmap 
#' 
#' This function is internal to 'pheatmap' and not exported, but 
#' our code depends on it. We added it here till it becomes exported.
#' 
#' @param x numeric vector
#' @param n number of breaks
#' @param center logical defaults to False

pheatmap_generate_breaks <- function (x, n, center = F) 
{
  if (center) {
    m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    res = seq(-m, m, length.out = n + 1)
  }
  else {
    res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 
                1)
  }
  return(res)
}