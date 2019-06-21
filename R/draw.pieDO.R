#' Draw a pie chart at specific xy coordinates with a specified number of slices and color for each slice.
#' 
#' Draw a pie graph evenly divided into polygons
#' each colored according to cols
#' the center of the pie chart is at x y
#' and the radius is rad
#' 
draw.pieDO <- function(x = 0.5, y = 0.5, radius = 1, cols = c("red", "green"), 
                       border.col = "black", node.border.lwd = 1, labels = NULL, 
                       edges = 200, label.cex = 1, xlim = NULL, ylim = NULL, 
                       add = TRUE){
  
  pies <- rep(1, length(cols))
  init.angle = 0
  
  if(is.null(labels)){labels <- rep("", length(cols))}
  
  pies <- c(0, cumsum(pies)/sum(pies))
  d.pie <- diff(pies)
  n.pie <- length(d.pie)
  
  
  twopi <- 2 * pi
  t2xy <- function(x,y,t){
    t2p <- twopi * t + init.angle * pi/180
    list(x = (radius * cos(t2p))+x, y = (radius * sin(t2p))+y)
  }
  
  
  if(!add){
    plot.new()
    if(is.null(xlim)){xlim = c(x-radius, x+radius)}
    if(is.null(ylim)){ylim = c(y-radius, y+radius)}
    plot.window(xlim = xlim, ylim = ylim)
  }
  for(i in 1:n.pie){
    n <- max(2, floor(edges * d.pie[i]))
    P <- t2xy(x = x, y = y, seq.int(pies[i], pies[i + 1], length.out = n))
    polygon(c(P$x, x), c(P$y, y), col = cols[i], border = border.col, lwd = node.border.lwd)
    text.x <- mean(c(min(P$x), max(P$x)))
    text.y <- mean(c(min(P$y), max(P$y)))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)){
      text(text.x, text.y, labels[i], xpd = TRUE, adj = 0.5, cex = label.cex)
    }
  }
  
}
