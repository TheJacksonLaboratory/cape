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
#' @param data_obj a \code{\link{Cape}} object
#' @param orientation string, ("vertical", "horizontal")
#' @param pos_col The color to use for positive main effects and interactions
#' must be one of "green", "purple", "red", "orange", "blue", "brown", "yellow", "gray"
#' see \code{\link{get_color}}
#' @param neg_col The color to use for negative main effects and interactions
#' takes the same values as pos_col.
#' @param light_dark Indicates whether pos_col, neg_col, and bg_col should be selected
#' from light colors ("l"), dark colors ("d") or the full spectrum from light to dark ("f")
#' @param pheno_labels Vector of phenotype names if other than what is stored in the
#' data object
#' @param cex_barplot_axis Size of axis for the bar plot
#' @param cex_barplot_labels Size of labels for the bar plot
#' @param cex_barplot_title Size of the barplot title
#' @param main Title for the plot. Defaults to "Eigentrait Contributions to Phenotypes"
#' @param cex_main Size of the overall title
#' @param main_x x shift for the overall title
#' @param main_y y shift for the overall title
#' @param cex_ET Size of the eigentrait labels
#' @param ET_label_x x shift for the eigentrait labels
#' @param ET_label_y y shift for the eigentrait labels
#' @param pheno_label_pos x shift for the trait labels
#' @param cex_pheno size of the trait labels
#' @param pheno_srt Rotation factor for the trait labels
#' @param percent_total_variance_x x shift for the percent total variance labels
#' @param percent_total_variance_y y shift for the percent total variance labels
#' @param cex_color_scale label size for the color scal 
#' @param cex_var_accounted size for the variance accounted for labels
#' @param var_accounted_x x shift for the variance accounted axis label
#' @param var_accounted_y x shift for the variance accounted axis label
#' @param show_var_accounted logical
#' @param just_selected_et logical
#'
#' @return \code{list("data_obj" = data_obj, "geno_obj" = geno_obj)}
#' 
#' @importFrom graphics barplot image
#' 
#' @export
plot_svd <- function(data_obj, orientation = c("vertical", "horizontal"), neg_col = "blue", 
                    pos_col =  "brown", light_dark = "f", pheno_labels = NULL, cex_barplot_axis = 1.7, 
                    cex_barplot_labels = 2, cex_barplot_title = 1.7, 
                    main = "Eigentrait Contributions to Phenotypes", 
                    cex_main = 2, main_x = 0.5, main_y = 0.5, cex_ET = 1.7, 
                    ET_label_x = 0.5, ET_label_y = 0.5, pheno_label_pos = 0.5, 
                    cex_pheno = 1.7, pheno_srt = 90, percent_total_variance_x = 0.5, 
                    percent_total_variance_y = 0.5, cex_color_scale = 1, cex_var_accounted = 2, 
                    var_accounted_x = 0, var_accounted_y = 0, show_var_accounted = FALSE,
                    just_selected_et = FALSE){
  
  oldPar <- par(no.readonly = TRUE)
  on.exit(oldPar)

  #test to see if there is a v in the orientation vector
  orient_test <- grep("v", orientation) 
  
  if(is.null(pheno_labels)){
    pheno_labels <- colnames(data_obj$pheno)
  }
  if(length(pheno_labels) != ncol(data_obj$pheno)){
    stop("pheno_labels needs to have as many elements as there are phenotypes.\n")
  }
  
  et <- data_obj$ET
  if(just_selected_et){
    et_which <- as.numeric(gsub("ET", "", colnames(et)))
  }else{
    et_which <- 1:length(pheno_labels)
  }
  
  #The contribution of each eigentrait to each phenotype
  #is encoded in the right singular vectors (eigenvectors)
  eigen_weights <- t(data_obj$right_singular_vectors)
  colnames(eigen_weights) <- colnames(data_obj$pheno)
  rownames(eigen_weights) <- paste0("ET", 1:nrow(eigen_weights))
  
  #use the singular values to calculate the variance
  #captured by each mode
  singular_vals <- data_obj$singular_values
  
  if(is.null(singular_vals)){
    stop("get_eigentraits() must be run before plot_svd()")
  }
  
  eigen_vals <- singular_vals^2
  var_accounted <- (eigen_vals/sum(eigen_vals))[et_which]
  
  
  #if the orientation is vertical, use the matrix as is.
  #otherwise, rotate the matrix, so the ETs are plotted 
  #in each row
  if(length(orient_test) > 0){
    rotated_eigen <- eigen_weights[et_which,]
  }else{
    rotated_eigen <- rotate_mat(eigen_weights[et_which,])
  }
  
  
  min_weight <- min(eigen_weights)
  max_weight <- max(eigen_weights)
  
  pos_cols <- get_color(pos_col, light_dark)
  neg_cols <- get_color(neg_col, light_dark)[3:1]
  
  mypal_pos <- colorRampPalette(pos_cols)
  mypal_neg <- colorRampPalette(neg_cols)
  
  
  color_levels <- seq(min_weight, max_weight, length=256)
  color_ramp <- c(mypal_neg(length(which(color_levels < 0))), mypal_pos(length(which(color_levels >= 0))))
    
  
  #plot the horizontal configuration
  if(length(orient_test) == 0){
    if(main != ""){
      layout_mat <- matrix(c(0,1,1,0,3,2,5,7,0,4,6,0), nrow = 3, ncol = 4, byrow = TRUE)
      layout(layout_mat, heights = c(0.15, 0.7, 0.15), widths = c(0.1, 0.35, 0.45, 0.1))
      #layout.show(7)
      
      #1) plot the title in its own box
      par(mar = c(0,0,0,0))
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))
      text(main_x, main_y, main, cex = cex_main)
    }else{
      layout_mat <- matrix(c(2,1,4,6,0,3,5,0), nrow = 2, ncol = 4, byrow = TRUE)
      layout(layout_mat, heights = c(0.7, 0.15), widths = c(0.1, 0.35, 0.45, 0.1))
      #layout.show(6)
      
    }
    
    #2) plot the weight matrix
    par(mar = c(0, 0, 2, 0))
    # my_image_plot(rotated_eigen)
    image(x = 1:length(rotated_eigen[,1]), y = 1:length(rotated_eigen[1,]), 
    rotated_eigen, col = color_ramp, axes = FALSE, xlab = "", ylab = "")
    
    #3) fit in the text for the y axis of the weight matrix
    par(mar = c(0,0,0,0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, ncol(rotated_eigen)))
    text(x = ET_label_x, y = seq(ET_label_y, (ncol(rotated_eigen)-0.5), 1), 
    rev(rownames(eigen_weights)[et_which]), cex = cex_ET) #the labels for the y axis (phenotypes)
    
    #4) fit in the text for the x axis of the weight matrix
    par(mar = c(0,0,0,0))
    plot.new()
    plot.window(ylim = c(0, 1), xlim = c(0, length(eigen_weights[1,])))
    text(y = pheno_label_pos, x = seq(0.5, (length(eigen_weights[1,])-0.5), 1), pheno_labels, cex = cex_pheno, srt = pheno_srt) #the labels for the y axis (phenotypes)
    
    
    #5) plot the barplot
    if(show_var_accounted){par(mar = c(0, 0, 2, 6))}else{par(mar = c(0,0,2,2))}
    a <- barplot(rev(var_accounted), names = "", horiz = TRUE, cex.axis = cex_barplot_axis, cex.lab = cex_barplot_labels, xlab = "")
    if(show_var_accounted){
      par(xpd = TRUE)
      text(y = a[,1], x = rev(var_accounted)+(max(rev(var_accounted))*0.01)+var_accounted_x, labels = paste0(signif(rev(var_accounted)*100, 2), "%"), adj = 0, cex = cex_var_accounted)
      par(xpd = FALSE)
    }
    
    #6) add the axis label for the barplot
    par(mar = c(0,0,0,0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    par(xpd = TRUE)
    text(x = percent_total_variance_x, y = percent_total_variance_y, labels = "% Total Var.", cex = cex_barplot_title) #the labels for the y axis (phenotypes)
    par(xpd = FALSE)		
    
    #7) add the color ramp
    # par(mar = c(0,2,2,2))
    par(mar = c(0,0,2,2))
    image(1, color_levels, matrix(data=color_levels, ncol=length(color_levels),nrow=1), col=color_ramp, xlab="",ylab="",xaxt="n", cex.axis = cex_color_scale)		
    
    # dev.off()
  }
  
  #If the orientation is set to vertical, there are lots of little parameters to adjust...
  if(length(orient_test) > 0){
    
    if(main != ""){
      layout_mat <- matrix(c(0,1,0,3,2,0,5,4,7,0,6,0), nrow = 4, ncol = 3, byrow = TRUE)
      layout(layout_mat, widths = c(0.2, 0.7, 0.1), heights = c(0.1, 0.35, 0.45, 0.1))
      #layout.show(7)
      
      #1) plot the title in its own box
      par(mar = c(0,0,0,0))
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))
      text(main_x, main_y, main, cex = cex_main)
    }else{
      layout_mat <- matrix(c(2,1,0,4,3,6,0,5,0), nrow = 3, ncol = 3, byrow = TRUE)
      layout(layout_mat, widths = c(0.2, 0.7, 0.1), heights = c(0.35, 0.45, 0.1))
      #layout.show(6)				
    }
    
    #2) plot the barplot
    par(mar = c(0, 0, 4, 2))
    a <- barplot(var_accounted, names = "", horiz = FALSE, cex.axis = cex_barplot_axis, cex.lab = cex_barplot_labels, ylab = "% Total Var.")
    if(show_var_accounted){
      par(xpd = TRUE)
      text(x = a[,1], y = var_accounted+(max(var_accounted)*0.05)+var_accounted_y, labels = paste0(signif(var_accounted*100, 2), "%"), adj = 0.5, cex = cex_var_accounted)
      par(xpd = FALSE)
    }
    
    
    #3) add the axis label for the barplot
    par(mar = c(0,0,0,0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    par(xpd = TRUE)
    text(x = percent_total_variance_x, y = percent_total_variance_y, 
    labels = "% Total Var.", cex = cex_barplot_title, srt = 90) #the labels for the y axis (phenotypes)
    par(xpd = FALSE)
    
    #3) plot the weight matrix
    par(mar = c(0, 0, 0, 2))
    image(x = 1:nrow(rotated_eigen), y = 1:ncol(rotated_eigen), rotated_eigen, 
    col = color_ramp, axes = FALSE, xlab = "", ylab = "")
    
    #4) fit in the text for the y axis of the weight matrix
    par(mar = c(0,0,0,0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, length(eigen_weights[1,])))
    text(x = pheno_label_pos, y = seq(0.5, (ncol(rotated_eigen)-0.5), 1), 
    labels = pheno_labels, cex = cex_pheno) #the labels for the y axis (phenotypes)
    
    #5) fit in the text for the x axis of the weight matrix
    par(mar = c(0,0,0,0))
    plot.new()
    plot.window(xlim = c(0, nrow(rotated_eigen)), ylim = c(0, 1))
    text(x = seq(ET_label_x, (nrow(rotated_eigen)-0.5), 1), 
    y = ET_label_y, rownames(rotated_eigen), cex = cex_ET) #the labels for the x axis (ETs)
    
    #6) add the color ramp
    #par(mar = c(0,2,0,2))
    par(mar = c(0,1,0,1))
    image(1, color_levels, matrix(data=color_levels, ncol=length(color_levels),nrow=1), col=color_ramp, xlab="",ylab="",xaxt="n", cex.axis = cex_color_scale)		
    
    names(var_accounted) <- colnames(data_obj$ET)
    invisible(var_accounted)
    
    # dev.off()			
    
  }
  
  
}
