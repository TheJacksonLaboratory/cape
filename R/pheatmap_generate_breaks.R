#' pheatmap generate breaks found at this link
#' https://cran.r-project.org/package=pheatmap 
#' 
#' This function is internal to 'pheatmap' and not exported, but 
#' our code depends on it. We added it here till it becomes exported.
#' 
#' @param x numeric vector
#' @param n number of breaks
#' @param center logical defaults to False
#' @keywords internal

pheatmap_generate_breaks <- function (x, n, center = FALSE) 
{
  if (center) {
    m = max(abs(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))))
    res = seq(-m, m, length.out = n + 1)
  }
  else {
    res = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 
                1)
  }
  return(res)
}