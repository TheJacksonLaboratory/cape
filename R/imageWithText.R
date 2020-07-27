#' Plot a heatmap
#'
#' @param mat A numeric matrix to plot as a heat map
#' @param xlab A string label for the x axis
#' @param ylab A string label for the y axis
#' @param main A title for the plot
#' @param main.shift A numeric value to shift the title
#' along the y axis.
#' @param col.names A vector of strings indicating names
#' for the columns of the matrix. Defaults to existing column
#' names.
#' @param row.names A vector of strings indicating names
#' for the rows of the matrix. Defaults to existing row
#' names.
#' @param row.text.adj adjustment value for row text indicating
#' the centering. See \link{\code{text}}.
#' @param row.text.shift numeric value for shifting the row labels
#' toward or away from the matrix.
#' @param row.text.rotation Rotation value in degrees for the row labels
#' @param col.text.rotation Rotation value in degrees for the column labels
#' @param col.text.adj adjustment value for column text indicating
#' the centering. See \link{\code{text}}.
#' @param col.text.shift numeric value for shifting the column labels
#' toward or away from the matrix.
#' @param show.text Whether to write the numerical value of each cell 
#' in the plot.
#' @param cex The size of the text when show.text is TRUE
#' @param col.text.cex The size of the column labels
#' @param row.text.cex The size of the row labels
#' @param main.cex The size of the title of the plot
#' @param split.at.vals Whether to split the values into
#' different color classes
#' @param split.points If split.at.vals is TRUE, split.points
#' determines the boundaries of the classes. For example, if
#' split.points is 0, negative and positive numbers will be 
#' assigned to different color classes.
#' @param col.scale One of c("green", "purple", "orange", "blue", 
#' "brown", "gray") to indicate the color scale to be used. Defaults
#' to gray.
#' @param color.spread A numerical value used as input to 
#' \link{\code{get.color2}} indicating the numeric distance
#' between colors in a ramp. Smaller values produce a smaller
#' difference between adjacent colors in the ramp.
#' @param light.dark One of "l", "d", or "f" indicating whether
#' the colors used should be light ("l"), dark ("d"), or from
#' across the full spectrum ("f").
#' @param class.mat An optional numeric matrix defining the color 
#' classes of each cell in the matrix. If omitted this is calculated
#' by the function.
#' @param grad.dir A string specifying how the color gradient 
#' should be applied. If "high" higher values are given darker
#' colors. If "low", lower values are given darker colors.
#' If "middle" values in the middle of the spectrum are 
#' given darker colors, and if "end" values at the ends
#' of the spectrum are given darker colors.
#' @param color.fun Either "linear" or "exponential" indicating
#' how the colors should transition from light to dark across
#' values.
#' @param exp.steepness If color.fun is "exponential," exp.steepness
#' indicates how quickly the colors should transition from light 
#' to dark.
#' @param global.color.scale Whether to impose a global minimum
#' and maximum to the colors, or to use the values themselves to 
#' determine the top and bottom of the color scale.
#' @param global.min If global.color.scale is TRUE, the minimum
#' value that should be assigned a color.
#' @param global.max If global.color.scale is TRUE, the maximum
#' value that should be assigned a color.
#' @param sig.digs The number of significant figures to use
#' from the input matrix. Helpful primarily when show.text is TRUE.
#' @param use.pheatmap.colors If TRUE, all other color parameters
#' are ignored, and colors like those used in the R package
#' pheatmap are used instead.
#' @param na.col The color to use for missing values.
#' @param gridlines Whether to plot gridlines on the matrix.
#'
#' @return None


col.names = colnames(mat), row.names = rownames(mat), row.text.adj = 1, row.text.shift = 0, 
row.text.rotation = 0, col.text.rotation = 90, col.text.adj = 1, 
col.text.shift = 0, show.text = TRUE, cex = 0.5, col.text.cex = 1, 
row.text.cex = 1, main.cex = 1, split.at.vals = FALSE, split.points = 0, 
col.scale = "gray", color.spread = 50, light.dark = "f", class.mat = NULL, 
grad.dir = c("high", "low", "middle", "ends"), color.fun = c("linear", "exponential"), 
exp.steepness = 1, global.color.scale = FALSE, global.min = NULL, global.max = NULL, 
sig.digs = 3, use.pheatmap.colors = FALSE, na.col = "lightgray", gridlines = FALSE

#' @param split.at.vals A logical value indicating whether
#' the numbers should be broken into classes with different 
#' colors.
#' @param split.points If split.at.vals is TRUE, split.points
#' is used to define multiple classes of numbers. For example,
#' if split.points is 0, negative numbers will be assigned one
#' class of colors, and positive numbers will be assigned another.
#' @param col.scale One of c("green", "purple", "orange", "blue", 
#' "brown", "gray") to indicate the color scale to be used. Defaults
#' to gray.
#' @param light.dark One of "l", "d", or "f" indicating whether
#' the colors used should be light ("l"), dark ("d"), or from
#' across the full spectrum ("f").
#' @param grad.dir A string specifying how the color gradient 
#' should be applied. If "high" higher values are given darker
#' colors. If "low", lower values are given darker colors.
#' If "middle" values in the middle of the spectrum are 
#' given darker colors, and if "end" values at the ends
#' of the spectrum are given darker colors.
#' @param color.fun Either "linear" or "exponential" indicating
#' how the colors should transition from light to dark across
#' values.
#' @param exp.steepness If color.fun is "exponential," exp.steepness
#' indicates how quickly the colors should transition from light 
#' to dark.
#' @param global.color.scale Whether to impose a global minimum
#' and maximum to the colors, or to use the values themselves to 
#' determine the top and bottom of the color scale.
#' @param global.min If global.color.scale is TRUE, the minimum
#' value that should be assigned a color.
#' @param global.max If global.color.scale is TRUE, the maximum
#' value that should be assigned a color.
#' @param use.pheatmap.colors If TRUE, all other color parameters
#' are ignored, and colors like those used in the R package
#' pheatmap are used instead.
#' @param na.col The color to use for missing values.
#'


imageWithText <- function(mat, xlab = "", ylab = "", main = NULL, main.shift = 0.12, 
col.names = colnames(mat), row.names = rownames(mat), row.text.adj = 1, row.text.shift = 0, 
row.text.rotation = 0, col.text.rotation = 90, col.text.adj = 1, 
col.text.shift = 0, show.text = TRUE, cex = 0.5, col.text.cex = 1, 
row.text.cex = 1, main.cex = 1, split.at.vals = FALSE, split.points = 0, 
col.scale = "gray", color.spread = 50, light.dark = "f", class.mat = NULL, 
grad.dir = c("high", "low", "middle", "ends"), color.fun = c("linear", "exponential"), 
exp.steepness = 1, global.color.scale = FALSE, global.min = NULL, global.max = NULL, 
sig.digs = 3, use.pheatmap.colors = FALSE, na.col = "lightgray", gridlines = FALSE){

		require(grid)
		
		#make sure Inf and -Inf are coded as NA
		mat[which(!is.finite(mat))] <- NA
		
		if(length(which(is.na(mat))) == length(mat)){
			return()
			}
		# if(length(light.dark) < length(col.scale)){light.dark <- rep(light.dark, length(col.scale))}
		
		get.default.col.fun <- grep("lin", color.fun)
		if(length(get.default.col.fun) > 0){
			color.fun = "linear"
			}
		
		end.fudge.factor = 10^-10

		if(is.null(class.mat)){
			class.mat <- matrix(1, dim(mat)[1], dim(mat)[2])
			}

		if(split.at.vals){
			for(p in 1:length(split.points)){
				class.mat[which(mat >= split.points[p])] <- class.mat[which(mat >= split.points[p])] + 1
				}
			# if(length(grad.dir) == 2){grad.dir <- "ends"}else{grad.dir <- "high"}		
			}else{
			split.points <- NULL	
			}
		class.mat[which(is.na(mat))] <- NA

		while(min(class.mat, na.rm = TRUE) > 1){class.mat <- class.mat - 1}
		
		classes <- sort(unique(as.vector(class.mat)))
		num.classes <- length(classes)
		if(num.classes == 1){
			class.mat <- matrix(1, nrow(mat), ncol(mat))
			}
	
			
		if(length(col.scale) == (length(split.points)+1)){
			class.cols <- col.scale
			}else{
			class.cols <- col.scale[classes]
			}
		if(length(col.scale) < num.classes){
			extra.cols <- num.classes - length(col.scale)
			col.scale <- c(col.scale, rep(col.scale, ceiling(extra.cols/length(col.scale))))
			}
			
		get.default <- grep("h", grad.dir)
		if(length(get.default) > 0){
			grad.dir <- "high"
			}
		
		
		# if(light.dark == "f"){max.col = 4}else{max.col = 8}
		max.col = 4
		dir.list <- vector(mode = "list", length = num.classes)
		names(dir.list) <- classes
		if(grad.dir == "high"){
			for(i in 1:length(dir.list)){
				dir.list[[i]] <- 1:max.col
				}
			}
		if(grad.dir == "low"){
			for(i in 1:length(dir.list)){
				dir.list[[i]] <- max.col:1
				}
			}
		if(grad.dir == "middle"){
			if(length(dir.list) != 2){stop("I can only color the middle if there are exactly two classes")}
			dir.list[[1]] <- 1:max.col
			dir.list[[2]] <- max.col:1
			}
			
			
		if(grad.dir == "ends"){
			if(length(dir.list) != 2){stop("I can only color the ends if there are exactly two classes")}
			dir.list[[1]] <- max.col:1
			dir.list[[2]] <- 1:max.col
			}



		#============================================================================
		#internal functions
		#============================================================================
		#This function takes in a matrix of values matched with colors, and 
		#a vector of values. It matched up the appropriate color for each
		#value in the vector
		bin.cols <- function(color.key, V){
			color.v <- rep(NA, length(V))
			for(i in 1:length(V)){
				diff.v  <- V[i] - as.numeric(color.key[,1])
				closest.val <- which(abs(diff.v) == min(abs(diff.v)))[1]
				color.v[i] <- color.key[closest.val,2]
				}
			return(color.v)
			}

		#This function generates the matrix of colors to use in 
		#the raster function
		fill.color.ramp <- function(mat, class.mat, global){

			#make color scales for each class or globally as defined
			color.scales <- vector(mode = "list", length = num.classes)
			names(color.scales) <- classes

			ColorRamp <- matrix(NA, dim(mat)[1], dim(mat)[2])
			num.classes = length(unique(as.vector(class.mat)))
			
			for(cl in 1:length(classes)){
				if(global.color.scale){
					if(is.null(global.min)){min.cl = min(mat, na.rm = TRUE)}else{min.cl = global.min}
					if(is.null(global.max)){max.cl = max(mat, na.rm = TRUE)}else{max.cl = global.max}	
					}else{
					if(length(which(class.mat == classes[cl])) > 0 && !all(is.na(mat[which(class.mat == classes[cl])]))){
						min.cl <- min(mat[which(class.mat == classes[cl])], na.rm = TRUE)
						max.cl <- max(mat[which(class.mat == classes[cl])], na.rm = TRUE)
						}else{
						min.cl <- NA	
						}
						
					}
					
				if(!is.na(min.cl)){
					if(color.fun == "linear"){
						ColorLevels <- seq(min.cl, max.cl, length=256)
						}else{
						ColorLevels <- exp.color.fun(min.cl, max.cl, steepness = exp.steepness, num.cols=256)	
						}
	
					#make the function to generate 
					col.vals <- get.color2(col.scale[cl], col.gap = color.spread)
										
					color.locale <- which(names(color.scales) == classes[cl])
					color.scales[[color.locale]] <- colorRampPalette(col.vals[dir.list[[color.locale]]])
					
					#find the entries in each class
					entry.locale <- which(class.mat == classes[cl])
					if(length(entry.locale) > 0){
						entry.vals <- mat[entry.locale]
						color.key <- cbind(ColorLevels, do.call(color.scales[[color.locale]], list(256)))
						entry.cols <- bin.cols(color.key, entry.vals)
						entry.order <- order(entry.locale)
						ColorRamp[entry.locale[entry.order]] <- entry.cols
						}
					}
				}
			
			return(ColorRamp)
			
			}
		
		#============================================================================


		if(use.pheatmap.colors){
			pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
			if(global.color.scale){
				bks <- pheatmap:::generate_breaks(c(global.min, global.max), length(pal), 
				center = F)
			}else{
				bks <- pheatmap:::generate_breaks(mat, length(pal), center = F)
			}
			ColorRamp <- pheatmap:::scale_colours(mat, col=pal, breaks=bks, na_col = na.col)
		}else{
			ColorRamp = fill.color.ramp(mat, class.mat, global.color.scale)
		}
		 	
		zmin <- min(mat, na.rm = TRUE); zmax <- max(mat, na.rm = TRUE)

		na.locale <- which(is.na(mat))
		if(length(na.locale) > 0){
			mat[na.locale] <- 0
			}
		
		all.ind <- which(!is.na(mat), arr.ind = TRUE)
		if(length(na.locale) > 0){
			mat[na.locale] <- NA
			}
		
		#translate the ColorRamp matrix to rgb matrices so we can use grid.
		rgb.mat <- col2rgb(ColorRamp)
		col <- matrix(rgb(rgb.mat["red",]/256,rgb.mat["green",]/256,rgb.mat["blue",]/256), dim(mat)[1], dim(mat)[2])

		max.dim <- max(c(dim(mat)[1], dim(mat)[2]))

		plot(c(1, dim(mat)[2]), c(1, dim(mat)[1]), type = "n", axes = FALSE, xlab = xlab, ylab = ylab, xlim = c(0.7, dim(mat)[2]+0.2), ylim = c(0.7, dim(mat)[1]+0.2), bg = "transparent")

		rasterImage(col, xleft = 0.5, ybottom = 0.5, xright = dim(mat)[2]+0.5, ytop = dim(mat)[1]+0.5, interpolate = FALSE, bg = "transparent")
		x.coord <- matrix(segment.region(0.5, dim(mat)[2]+0.5, dim(mat)[2], "center"), nrow = dim(mat)[1], ncol = dim(mat)[2], byrow = TRUE)
		y.coord <- matrix(segment.region(0.5, dim(mat)[1]+0.5, dim(mat)[1], "center"), nrow = dim(mat)[1], ncol = dim(mat)[2])
		
		if(show.text){
		text(x.coord, rev(y.coord), labels = signif(as.vector(mat), sig.digs), cex = cex)
		}
		
		if(!is.null(main)){
			par(xpd = TRUE)
			plot.range = max(y.coord) - min(y.coord)
			text(mean(x.coord[1,]), max(y.coord)+(plot.range*main.shift), labels = main, cex = main.cex)
			par(xpd = FALSE)
			}
			
		if(!is.null(col.names)){
			if(length(col.names) != dim(mat)[2]){
				stop("There is a different number of column names than columns.")
				}
			par(xpd = TRUE)
			# text(x.coord[1,], (min(y.coord)-(max(y.coord)) - col.text.shift), labels = col.names, srt = col.text.rotation, adj = col.text.adj)
			plot.range <- max(y.coord) - min(y.coord)
			text((x.coord[1,]), (min(y.coord) - (plot.range*0.1) - col.text.shift), labels = col.names, srt = col.text.rotation, adj = col.text.adj, cex = col.text.cex)
			}

		if(!is.null(row.names)){
			if(length(row.names) != dim(mat)[1]){
				stop("There is a different number of row names than rows.")
				}
			par(xpd = TRUE)
			plot.range <- (max(x.coord) - min(x.coord))
			text((min(x.coord) - (plot.range*0.1) - row.text.shift), y.coord[,1], labels = rev(row.names), adj = row.text.adj, cex = row.text.cex, srt = row.text.rotation)
			}


	}