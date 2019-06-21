#' This function plots the trait circles which have gaps at
#' the right-hand sice with labels
plot.trait.circ <- function(trait.circ, label.gap, plot.dim, main.lwd){
  
  num.pheno = length(trait.circ)
  
  pheno.label.starts <- round(segment.region(1, label.gap, num.points = num.pheno, alignment = "center"))
  
  #put the labels between the outermost trait circle and the edge of the plot
  label.x <- max(trait.circ[[length(trait.circ)]]$x) + (plot.dim[2] - max(trait.circ[[length(trait.circ)]]$x))*0.25
  arrow.start.x <- max(trait.circ[[length(trait.circ)]]$x) + (plot.dim[2] - max(trait.circ[[length(trait.circ)]]$x))*0.2
  
  for(tr in length(trait.circ):1){
    #add light gray bars to help the eye track the phenotypes
    #they should be staggered so the label sticks don't overlap
    circ.x <- trait.circ[[tr]]$x[pheno.label.starts[tr]:length(trait.circ[[tr]]$x)]
    circ.y <- trait.circ[[tr]]$y[pheno.label.starts[tr]:length(trait.circ[[tr]]$y)]
    points(circ.x, circ.y, type = "l", col = "lightgray", lwd = main.lwd)
    
    #and add phenotype labels
    arrow.end.x <- trait.circ[[tr]]$x[pheno.label.starts[tr]]
    arrow.y <- trait.circ[[tr]]$y[pheno.label.starts[tr]]
    segments(x0 = arrow.start.x, x1 = arrow.end.x, y0 = arrow.y, lwd = main.lwd, col = "lightgray")
    text(label.x, arrow.y, labels = names(trait.circ)[tr], adj = 0)
  }
}
