#' Plot concentric trait circles
#'
#' This internal function is called by \code{\link{plot_network}}
#' It plots one circle for each trait in the data object.
#'
#' @param trait_circ A list with one element for each cirlce giving
#' the coordinates of the points to plot for each circle.
#' @param label_gap The amount of space to leave between the circles
#' and their labels
#' @param plot_dim the dimensions of the plot
#' @param main_lwd The line thickness for the circles
#'
#' @return None
#' @keywords internal

plot_trait_circ <- function(trait_circ, label_gap, plot_dim, main_lwd){
  
  num_pheno = length(trait_circ)
  
  pheno_label_starts <- round(segment_region(1, label_gap, num_points = num_pheno, alignment = "center"))
  
  #put the labels between the outermost trait circle and the edge of the plot
  label_x <- max(trait_circ[[length(trait_circ)]]$x) + (plot_dim[2] - max(trait_circ[[length(trait_circ)]]$x))*0.25
  arrow_start_x <- max(trait_circ[[length(trait_circ)]]$x) + (plot_dim[2] - max(trait_circ[[length(trait_circ)]]$x))*0.2
  
  for(tr in length(trait_circ):1){
    #add light gray bars to help the eye track the phenotypes
    #they should be staggered so the label sticks don't overlap
    circ_x <- trait_circ[[tr]]$x[pheno_label_starts[tr]:length(trait_circ[[tr]]$x)]
    circ_y <- trait_circ[[tr]]$y[pheno_label_starts[tr]:length(trait_circ[[tr]]$y)]
    points(circ_x, circ_y, type = "l", col = "lightgray", lwd = main_lwd)
    
    #and add phenotype labels
    arrow_end_x <- trait_circ[[tr]]$x[pheno_label_starts[tr]]
    arrow_y <- trait_circ[[tr]]$y[pheno_label_starts[tr]]
    segments(x0 = arrow_start_x, x1 = arrow_end_x, y0 = arrow_y, lwd = main_lwd, col = "lightgray")
    text(label_x, arrow_y, labels = names(trait_circ)[tr], adj = 0)
  }
}
