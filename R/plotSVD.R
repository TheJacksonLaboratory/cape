#' Plots eigentraits
#' 
#' This function plots the results of the singular value
#' decomposition (SVD) on the phenotypes. Gray bars indicate
#' the amount of phenotypic variance accounted for by each
#' eigentrait.
#' 
#' Below the bars is a heatmap indicating how each trait 
#' contributes to each eigentrait. Colors can be adjusted
#' to suit preferences.
#'
#' @param data.obj a \code{\link{Cape}} object
#' @param orientation string, ("vertical", "horizontal")
#' @param pos.col The color to use for positive main effects and interactions
#' must be one of "green", "purple", "red", "orange", "blue", "brown", "yellow", "gray"
#' see \code{\link{get.color}}
#' @param neg.col The color to use for negative main effects and interactions
#' takes the same values as pos.col.
#' @param light.dark Indicates whether pos.col, neg.col, and bg.col should be selected
#' from light colors ("l"), dark colors ("d") or the full spectrum from light to dark ("f")
#' @param pheno.labels Vector of phenotype names if other than what is stored in the
#' data object
#' @param cex.barplot.axis Size of axis for the bar plot
#' @param cex.barplot.labels Size of labels for the bar plot
#' @param cex.barplot.title Size of the barplot title
#' @param main Title for the plot. Defaults to "Eigentrait Contributions to Phenotypes"
#' @param cex.main Size of the overall title
#' @param main.x x shift for the overall title
#' @param main.y y shift for the overall title
#' @param cex.ET Size of the eigentrait labels
#' @param ET.label.x x shift for the eigentrait labels
#' @param ET.label.y y shift for the eigentrait labels
#' @param pheno.label.pos x shift for the trait labels
#' @param cex.pheno size of the trait labels
#' @param pheno.srt Rotation factor for the trait labels
#' @param percent.total.variance.x x shift for the percent total variance labels
#' @param percent.total.variance.y y shift for the percent total variance labels
#' @param cex.color.scale label size for the color scal 
#' @param cex.var.accounted size for the variance accounted for labels
#' @param var.accounted.x x shift for the variance accounted axis label
#' @param var.accounted.y x shift for the variance accounted axis label
#' @param show.var.accounted logical
#'
#' @return \code{list("data.obj" = data.obj, "geno.obj" = geno.obj)}
#'
#' @export
plotSVD <- function(data.obj, orientation = c("vertical", "horizontal"), neg.col = "blue", 
                    pos.col =  "brown", light.dark = "f", pheno.labels = NULL, cex.barplot.axis = 1.7, 
                    cex.barplot.labels = 2, cex.barplot.title = 1.7, 
                    main = "Eigentrait Contributions to Phenotypes", 
                    cex.main = 2, main.x = 0.5, main.y = 0.5, cex.ET = 1.7, 
                    ET.label.x = 0.5, ET.label.y = 0.5, pheno.label.pos = 0.5, 
                    cex.pheno = 1.7, pheno.srt = 90, percent.total.variance.x = 0.5, 
                    percent.total.variance.y = 0.5, cex.color.scale = 1, cex.var.accounted = 2, 
                    var.accounted.x = 0, var.accounted.y = 0, show.var.accounted = FALSE
){
  
  #test to see if there is a v in the orientation vector
  orient.test <- grep("v", orientation) 
  
  if(is.null(pheno.labels)){
    pheno.labels <- colnames(data.obj$pheno)
  }
  if(length(pheno.labels) != ncol(data.obj$pheno)){
    stop("pheno.labels needs to have as many elements as there are phenotypes.\n")
  }
  
  
  #The contribution of each eigentrait to each phenotype
  #is encoded in the right singular vectors (eigenvectors)
  eigen.weights <- t(data.obj$right_singular_vectors)[1:dim(data.obj$ET)[2],]
  colnames(eigen.weights) <- colnames(data.obj$pheno)
  rownames(eigen.weights) <- colnames(data.obj$ET)
  
  #use the singular values to calculate the variance
  #captured by each mode
  singular.vals <- data.obj$singular_values[1:dim(data.obj$ET)[2]]
  
  if(is.null(singular.vals)){
    stop("get.eigentraits() must be run before plotSVD()")
  }
  
  eigen.vals <- singular.vals^2
  var.accounted <- eigen.vals/sum(eigen.vals)
  
  
  #if the orientation is vertical, use the matrix as is.
  #otherwise, rotate the matrix, so the ETs are plotted 
  #in each row
  if(length(orient.test) > 0){
    rotated.eigen <- eigen.weights
  }else{
    rotated.eigen <- rotate.mat(eigen.weights)
  }
  
  
  min.weight <- min(eigen.weights)
  max.weight <- max(eigen.weights)
  
  pos.cols <- get.color(pos.col, light.dark)
  neg.cols <- get.color(neg.col, light.dark)[3:1]
  
  mypal.pos <- colorRampPalette(pos.cols)
  mypal.neg <- colorRampPalette(neg.cols)
  
  
  ColorLevels <- seq(min.weight, max.weight, length=256)
  ColorRamp <- c(mypal.neg(length(which(ColorLevels < 0))), mypal.pos(length(which(ColorLevels >= 0))))
  
  
  
  #plot the horizontal configuration
  if(length(orient.test) == 0){
    if(main != ""){
      layout.mat <- matrix(c(0,1,1,0,3,2,5,7,0,4,6,0), nrow = 3, ncol = 4, byrow = TRUE)
      layout(layout.mat, heights = c(0.15, 0.7, 0.15), widths = c(0.1, 0.35, 0.45, 0.1))
      #layout.show(7)
      
      #1) plot the title in its own box
      par(mar = c(0,0,0,0))
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))
      text(main.x, main.y, main, cex = cex.main)
    }else{
      layout.mat <- matrix(c(2,1,4,6,0,3,5,0), nrow = 2, ncol = 4, byrow = TRUE)
      layout(layout.mat, heights = c(0.7, 0.15), widths = c(0.1, 0.35, 0.45, 0.1))
      #layout.show(6)
      
    }
    
    #2) plot the weight matrix
    par(mar = c(0, 0, 2, 0))
    # myImagePlot(rotated.eigen)
    image(x = 1:length(rotated.eigen[,1]), y = 1:length(rotated.eigen[1,]), rotated.eigen, col = ColorRamp, axes = FALSE, xlab = "", ylab = "")
    
    #3) fit in the text for the y axis of the weight matrix
    par(mar = c(0,0,0,0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, length(eigen.weights[,1])))
    text(x = ET.label.x, y = seq(ET.label.y, (length(eigen.weights[,1])-0.5), 1), rev(rownames(eigen.weights)), cex = cex.ET) #the labels for the y axis (phenotypes)
    
    #4) fit in the text for the x axis of the weight matrix
    par(mar = c(0,0,0,0))
    plot.new()
    plot.window(ylim = c(0, 1), xlim = c(0, length(eigen.weights[1,])))
    text(y = pheno.label.pos, x = seq(0.5, (length(eigen.weights[1,])-0.5), 1), pheno.labels, cex = cex.pheno, srt = pheno.srt) #the labels for the y axis (phenotypes)
    
    
    #5) plot the barplot
    if(show.var.accounted){par(mar = c(0, 0, 2, 6))}else{par(mar = c(0,0,2,2))}
    a <- barplot(rev(var.accounted), names = "", horiz = TRUE, cex.axis = cex.barplot.axis, cex.lab = cex.barplot.labels, xlab = "")
    if(show.var.accounted){
      par(xpd = TRUE)
      text(y = a[,1], x = rev(var.accounted)+(max(rev(var.accounted))*0.01)+var.accounted.x, labels = paste0(signif(rev(var.accounted)*100, 2), "%"), adj = 0, cex = cex.var.accounted)
      par(xpd = FALSE)
    }
    
    #6) add the axis label for the barplot
    par(mar = c(0,0,0,0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    par(xpd = TRUE)
    text(x = percent.total.variance.x, y = percent.total.variance.y, labels = "% Total Var.", cex = cex.barplot.title) #the labels for the y axis (phenotypes)
    par(xpd = FALSE)		
    
    #7) add the color ramp
    # par(mar = c(0,2,2,2))
    par(mar = c(0,0,2,2))
    image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp, xlab="",ylab="",xaxt="n", cex.axis = cex.color.scale)		
    
    # dev.off()
  }
  
  #If the orientation is set to vertical, there are lots of little parameters to adjust...
  if(length(orient.test) > 0){
    
    if(main != ""){
      layout.mat <- matrix(c(0,1,0,3,2,0,5,4,7,0,6,0), nrow = 4, ncol = 3, byrow = TRUE)
      layout(layout.mat, widths = c(0.2, 0.7, 0.1), heights = c(0.1, 0.35, 0.45, 0.1))
      #layout.show(7)
      
      #1) plot the title in its own box
      par(mar = c(0,0,0,0))
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))
      text(main.x, main.y, main, cex = cex.main)
    }else{
      layout.mat <- matrix(c(2,1,0,4,3,6,0,5,0), nrow = 3, ncol = 3, byrow = TRUE)
      layout(layout.mat, widths = c(0.2, 0.7, 0.1), heights = c(0.35, 0.45, 0.1))
      #layout.show(6)				
    }
    
    #2) plot the barplot
    par(mar = c(0, 0, 4, 2))
    a <- barplot(var.accounted, names = "", horiz = FALSE, cex.axis = cex.barplot.axis, cex.lab = cex.barplot.labels, ylab = "% Total Var.")
    if(show.var.accounted){
      par(xpd = TRUE)
      text(x = a[,1], y = var.accounted+(max(var.accounted)*0.05)+var.accounted.y, labels = paste0(signif(var.accounted*100, 2), "%"), adj = 0.5, cex = cex.var.accounted)
      par(xpd = FALSE)
    }
    
    
    #3) add the axis label for the barplot
    par(mar = c(0,0,0,0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    par(xpd = TRUE)
    text(x = percent.total.variance.x, y = percent.total.variance.y, labels = "% Total Var.", cex = cex.barplot.title, srt = 90) #the labels for the y axis (phenotypes)
    par(xpd = FALSE)
    
    #3) plot the weight matrix
    par(mar = c(0, 0, 0, 2))
    image(x = 1:length(eigen.weights[,1]), y = 1:length(eigen.weights[1,]), rotated.eigen, col = ColorRamp, axes = FALSE, xlab = "", ylab = "")
    
    #4) fit in the text for the y axis of the weight matrix
    par(mar = c(0,0,0,0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, length(eigen.weights[1,])))
    text(x = pheno.label.pos, y = seq(0.5, (length(eigen.weights[1,])-0.5), 1), labels = pheno.labels, cex = cex.pheno) #the labels for the y axis (phenotypes)
    
    #5) fit in the text for the x axis of the weight matrix
    par(mar = c(0,0,0,0))
    plot.new()
    plot.window(xlim = c(0, length(eigen.weights[,1])), ylim = c(0, 1))
    text(x = seq(ET.label.x, (length(eigen.weights[,1])-0.5), 1), y = ET.label.y, rownames(eigen.weights), cex = cex.ET) #the labels for the x axis (ETs)
    
    #6) add the color ramp
    #par(mar = c(0,2,0,2))
    par(mar = c(0,1,0,1))
    image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp, xlab="",ylab="",xaxt="n", cex.axis = cex.color.scale)		
    
    names(var.accounted) <- colnames(data.obj$ET)
    invisible(var.accounted)
    
    # dev.off()			
    
  }
  
  
}


