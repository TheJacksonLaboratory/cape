#' Generate a Heatmap-type image
#'
#' This internal function generates heatmap-type
#' images for functions like \code{\link{plot_svd}}
#' \code{\link{plot_pairscan}}, and \code{\link{plot_variant_influences}}
#' 
#' @param x Matrix to be plotted 

#' @param ... possible parameters are main, xlab, ylab, mark_coords,
#'  mark_col, show_labels, chromosome_coordinates, chr_names, chr_labels,
#'  show_pheno_labels, extra_col_mat,allele_cols
#'  
#' @importFrom stats median
#' @keywords internal

my_image_plot <- function(x,...){
    # print(dim(x))
    
    oldPar <- par(no.readonly = TRUE)
		on.exit(oldPar)

    #build the argument list from additional arguments added to
    #the function
    additional_arguments <- list(...)
    
    #=================================================================
    #There are some special additional arguments that this
    #function can take in
    #if there are xlab and ylab arguments here
    #pull them out first before we override them
    xlab <- additional_arguments$xlab
    ylab <- additional_arguments$ylab
    show_labels <- additional_arguments$show_labels
    additional_arguments$show_labels <- NULL
    if(is.null(show_labels)){show_labels <- TRUE}
    
    show_pheno_labels <- additional_arguments$show_pheno_labels
    additional_arguments$show_pheno_labels <- NULL
    if(is.null(show_pheno_labels)){show_pheno_labels <- TRUE}
    
    chromosome_coordinates <- additional_arguments$chromosome_coordinates
    additional_arguments$chromosome_coordinates <- NULL
    chr_names <- additional_arguments$chr_names
    additional_arguments$chr_names <- NULL
    chr_labels <- additional_arguments$chr_labels
    additional_arguments$chr_labels <- NULL
    
    mark_coords <- additional_arguments$mark_coords
    mark_col <- additional_arguments$mark_col
    additional_arguments$mark_coords <- NULL
    additional_arguments$mark_col <- NULL
    
    pos_col <- additional_arguments$pos_col
    if(is.null(pos_col)){pos_col = "brown"}
    additional_arguments$pos_col <- NULL
    neg_col <- additional_arguments$neg_col
    if(is.null(neg_col)){neg_col = "blue"}
    additional_arguments$neg_col <- NULL
    col_pal <- additional_arguments$col_pal
    if(is.null(col_pal)){col_pal = "f"}
    additional_arguments$col_pal <- NULL
    
    col_split_point <- additional_arguments$col_split_point
    if(is.null(col_split_point)){col_split_point <- 0}
    additional_arguments$col_split_point <- NULL
    
    allele_cols <- additional_arguments$allele_cols
    additional_arguments$allele_cols <- NULL
    
    extra_col_mat <- additional_arguments$extra_col_mat
    additional_arguments$extra_col_mat <- NULL
    
    #if there are min_x and max_x argments
    #pull these out too.
    min_x <- additional_arguments$min_x
    max_x <- additional_arguments$max_x
    #=================================================================
    
    
    #if they aren't specified, use the
    #x matrix to specify them
    if(is.null(min_x)){
      min_x <- min(x, na.rm = TRUE)
    }else{ #otherwise remove it from the argument list, so it doesn't throw a warning when we use it in image
      additional_arguments$min_x <- NULL
    }
    if(is.null(max_x)){
      max_x <- max(x, na.rm = TRUE)
    }else{
      additional_arguments$max_x <- NULL
    }
    
    
    yLabels <- rownames(x)
    xLabels <- colnames(x)
    
    layout_mat <- matrix(c(1:3, 4, 4, 4), nrow=2, ncol=3, byrow = TRUE)
    layout(layout_mat, widths=c(0.75,4,0.75), heights = c(1,0.1))
    
    color_levels <- seq(min_x, max_x, length=256)    
    
    pos_cols <- get_color(pos_col, col_pal)
    neg_cols <- get_color(neg_col, col_pal)[3:1]
    
    mypal_pos <- colorRampPalette(pos_cols)
    mypal_neg <- colorRampPalette(neg_cols)
    
    color_ramp <- c(mypal_neg(length(which(color_levels < col_split_point))), mypal_pos(length(which(color_levels >= col_split_point))))
    
    
    #=====================================
    #plot the y axis label
    par(mar = c(0,0,5,0))
    plot.new()
    plot.window(xlim = c(0,1), ylim = c(0,1))
    if(show_labels){
      text(x = 0.3, y = 0.5, ylab, srt = 90, cex = 2)
    }else{
      text(x = 0.85, y = 0.5, ylab, srt = 90, cex = 2)	
    }
    
    if(show_labels || show_pheno_labels){
      par(mar = c(5,3,5,2))
    }else{
      par(mar = c(3,3,5,2))
    }
    
    #add the default arguments to the argument list
    additional_arguments$x <- 1:length(xLabels)
    additional_arguments$y <- 1:length(yLabels)
    additional_arguments$z <- rotate_mat(x)
    # additional_arguments$z <- x
    additional_arguments$col = color_ramp
    additional_arguments$xlab <- ""
    additional_arguments$ylab = ""
    additional_arguments$axes = FALSE
    additional_arguments$zlim <- c(min_x, max_x)
    # additional_arguments$useRaster <- TRUE
    do.call(image, additional_arguments)
    
    #add the extra colors if we are highlighting particular cells
    if(!is.null(extra_col_mat)){
      u_col <- unique(as.vector(extra_col_mat[which(!is.na(extra_col_mat))]))
      if(length(u_col) > 0){
        for(i in 1:length(u_col)){
          col_locale <- which(rotate_mat(extra_col_mat) == u_col[i], arr.ind = TRUE)
          points(col_locale[,1], col_locale[,2], col = u_col[i], pch = 15, cex = 0.6)
        }
      }
    }
    
    chr_cols <- rep(c("darkgray", "white"), ceiling(length(chromosome_coordinates)/2))
    y_chr_coord <- length(yLabels) - chromosome_coordinates + 1
    poly_perc = 0.03; label_size <- 0.9
    plot_dim <- par("usr")
    poly_width <- plot_dim[2]*poly_perc
    poly_max <- plot_dim[1]; poly_min <- poly_max-poly_width; poly_mid <- mean(c(poly_min, poly_max))
    
    if(!is.null(chromosome_coordinates)){
      par(xpd = TRUE)
      for(i in 1:(length(chromosome_coordinates)-1)){
        if(dim(x)[2]+1 >= max(chromosome_coordinates)){
          polygon(x = c(chromosome_coordinates[i], chromosome_coordinates[i+1], chromosome_coordinates[i+1], chromosome_coordinates[i]), y = c(poly_min, poly_min, poly_max, poly_max), col = chr_cols[i])
          #if(chr_labels[i] == 0){chr_srt = 90}else{chr_srt = 0}
          chr_srt = 0
          #cat(chr_names[i], chr_labels[i], chr_srt, "\n")
          text(x = mean(c(chromosome_coordinates[i], chromosome_coordinates[i+1])), y = poly_mid, cex = label_size, labels = chr_names[i], srt = chr_srt)
        }
        chr_srt = 90
        #if(chr_labels[i] == 0){chr_srt = 0}else{chr_srt = 90}
        polygon(y = c(y_chr_coord[i], y_chr_coord[i+1], y_chr_coord[i+1], y_chr_coord[i]), x = c(poly_min, poly_min, poly_max, poly_max), col = chr_cols[i])
        text(y = mean(c(y_chr_coord[i], y_chr_coord[i+1])), x = poly_mid, cex = label_size, labels = chr_names[i], srt = chr_srt)
        
      }
      
      par(xpd = FALSE)
    }
    
    
    if(!is.null(allele_cols)){
      par(xpd = TRUE)
      for(i in 1:(length(chromosome_coordinates)-1)){
        text(x = 1:dim(x)[1], y = poly_min*(1.5+poly_perc), labels = rep("|", dim(x)[1]), col = allele_cols)
        # text(y = 1:dim(x)[1], x = poly_min*(1.2+poly_perc), labels = rep("|", dim(x)[1]), col = allele_cols)
      }
      par(xpd = FALSE)
    }
    
    
    if(!is.null(mark_coords)){
      points(mark_coords[,1], mark_coords[,2], col = mark_col, pch = 16, cex = 0.5)
    }
    
    
    if(show_labels){
      par(xpd = TRUE)
      if(is.null(chromosome_coordinates)){
        axis(LEFT<-2, at=1:length(yLabels), labels=rev(yLabels), las= HORIZONTAL<-1,cex.axis=0.7, adj = 1, tick = TRUE)
        axis(BELOW<-1, at=1:length(xLabels), labels=FALSE, cex.axis=0.7, tick = TRUE)
        text(x = 1:length(xLabels), y = poly_min*1.2, labels = xLabels, srt = 90, cex = 0.7, adj = 1)
      }else{
        axis(LEFT<-2, at=1:length(yLabels), labels=rev(yLabels), las= HORIZONTAL<-1,cex.axis=0.7, adj = 1, tick = FALSE, line = 0.8)
        axis(BELOW<-1, at=1:length(xLabels), labels=FALSE, cex.axis=0.7, tick = FALSE)
        text(x = 1:length(xLabels), y = poly_min*1.2, labels = xLabels, srt = 90, cex = 0.7, adj = 1)
      }
      par(xpd = FALSE)
    }else{
      axis(LEFT<-2, at=1:length(yLabels), labels=FALSE, las= HORIZONTAL<-1,cex.axis=0.7, adj = 1, tick = TRUE, lwd.ticks = 0)
      axis(BELOW<-1, at=1:length(xLabels), labels=FALSE, cex.axis=0.7, tick = TRUE, lwd.ticks = 0)
    }
    
    if(show_pheno_labels & !show_labels){
      only_pheno_labels <- xLabels
      marker_locale <- which(yLabels %in% xLabels)
      only_pheno_labels[marker_locale] <- ""
      par(xpd = TRUE)
      text(x = 1:length(only_pheno_labels), y = mean(c(poly_min, poly_max)), labels = only_pheno_labels, srt = 90, cex = 1.5, adj = 1)
      par(xpd = FALSE)
    }
    
    #plot the x axis label close to the axis if we are not printing the labels
    if(!show_labels){
      par(xpd = TRUE)
      plot_height = plot_dim[2]-plot_dim[1]
      marker_locale <- which(yLabels %in% xLabels)
      text(x = median(marker_locale), y = poly_min-(plot_height*0.05), xlab, cex = 2)
      par(xpd = FALSE)
    }
    
    
    
    # Color Scale
    par(mar = c(3,2.5,5,2))
    image(1, color_levels, matrix(data=color_levels, ncol=length(color_levels),nrow=1), col=color_ramp, xlab="",ylab="",xaxt="n", cex.axis = 2)
    
    #plot the x axis labels in a different window if we are printing the marker labels
    
    par(mar = c(0,0,0,2))
    plot.new()
    plot.window(xlim = c(0,1), ylim = c(0,1))
    if(show_labels){
      par(xpd = TRUE)
      text(x = 0.5, y = 0.5, xlab, cex = 2)
      par(xpd = FALSE)
    }
    
  }
