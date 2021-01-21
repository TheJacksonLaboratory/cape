#' Draw a pie chart
#' 
#' Draw a pie graph evenly divided into polygons
#' each colored according to the specified colors (cols)
#' The center of the pie chart is at coordinates x, y
#' and the radius is defined by rad
#' 
#' @param x the x coordinate for the center of the pie chart
#' @param y the y coordinate for the center of the pie chart
#' @param radius The radius of the pie chart
#' @param cols A vector of colors in any format specifying the colors 
#' of the sections of the pie chart. The number of colors determines
#' the number of sections in the pie.
#' @param border_col The color of the border of the pie chart. Use NA
#' to omit the border.
#' @param node_border_lwd The thickness of the border line
#' @param labels A vector of character values giving labels for the 
#' sections of the pie.
#' @param edges The default value is 200, which generates a circle.
#' Reducing this number will result in other shapes that may or may
#' not produce meaningful graphs depending on the number of pie sections.
#' @param label_cex A numeric value giving the size of the label text.
#' @param xlim A vector of length two giving the xlim of the plot if 
#' add is FALSE
#' @param ylim A vector of length two giving the ylim of the plot if 
#' add is FALSE
#' @param add A logical value indicating whether the pie should be 
#' added to an existing plot (TRUE) or drawn on a new plot (FALSE).
#' 
#' @importFrom graphics text
#' 
#' @keywords internal
#' 
draw_pie <- function(x = 0.5, y = 0.5, radius = 1, cols = c("red", "green"), border_col = "black", node_border_lwd = 1, 
labels = NULL, edges = 200, label_cex = 1, xlim = NULL, ylim = NULL, add = TRUE){
  
  pies <- rep(1, length(cols))
  init_angle = 0
  
  if(is.null(labels)){labels <- rep("", length(cols))}
  
  pies <- c(0, cumsum(pies)/sum(pies))
  d_pie <- diff(pies)
  n_pie <- length(d_pie)
  
  
  twopi <- 2 * pi
  t2xy <- function(x,y,t){
    t2p <- twopi * t + init_angle * pi/180
    list(x = (radius * cos(t2p))+x, y = (radius * sin(t2p))+y)
  }
  
  
  if(!add){
    plot.new()
    if(is.null(xlim)){xlim = c(x-radius, x+radius)}
    if(is.null(ylim)){ylim = c(y-radius, y+radius)}
    plot.window(xlim = xlim, ylim = ylim)
  }
  for(i in 1:n_pie){
    n <- max(2, floor(edges * d_pie[i]))
    P <- t2xy(x = x, y = y, seq.int(pies[i], pies[i + 1], length.out = n))
    polygon(c(P$x, x), c(P$y, y), col = cols[i], border = border_col, lwd = node_border_lwd)
    text_x <- mean(c(min(P$x), max(P$x)))
    text_y <- mean(c(min(P$y), max(P$y)))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)){
      text(text_x, text_y, labels[i], xpd = TRUE, adj = 0.5, cex = label_cex)
    }
  }
  
}
