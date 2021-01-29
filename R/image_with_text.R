#' Plot a heatmap
#'
#' @param mat A numeric matrix to plot as a heat map
#' @param xlab A string label for the x axis
#' @param ylab A string label for the y axis
#' @param main A title for the plot
#' @param main_shift A numeric value to shift the title
#' along the y axis.
#' @param col_names A vector of strings indicating names
#' for the columns of the matrix. Defaults to existing column
#' names.
#' @param row_names A vector of strings indicating names
#' for the rows of the matrix. Defaults to existing row
#' names.
#' @param row_text_adj adjustment value for row text indicating
#' the centering. See \code{\link{text}}.
#' @param row_text_shift numeric value for shifting the row labels
#' toward or away from the matrix.
#' @param row_text_rotation Rotation value in degrees for the row labels
#' @param col_text_rotation Rotation value in degrees for the column labels
#' @param col_text_adj adjustment value for column text indicating
#' the centering. See \code{\link{text}}.
#' @param col_text_shift numeric value for shifting the column labels
#' toward or away from the matrix.
#' @param show_text Whether to write the numerical value of each cell 
#' in the plot.
#' @param cex The size of the text when show_text is TRUE
#' @param col_text_cex The size of the column labels
#' @param row_text_cex The size of the row labels
#' @param main_cex The size of the title of the plot
#' @param split_at_vals Whether to split the values into
#' different color classes
#' @param split_points If split_at_vals is TRUE, split_points
#' determines the boundaries of the classes. For example, if
#' split_points is 0, negative and positive numbers will be 
#' assigned to different color classes.
#' @param col_scale One of c("green", "purple", "orange", "blue", 
#' "brown", "gray") to indicate the color scale to be used. Defaults
#' to gray.
#' @param color_spread A numerical value used as input to 
#' \code{\link{get_color2}} indicating the numeric distance
#' between colors in a ramp. Smaller values produce a smaller
#' difference between adjacent colors in the ramp.
#' @param light_dark One of "l", "d", or "f" indicating whether
#' the colors used should be light ("l"), dark ("d"), or from
#' across the full spectrum ("f").
#' @param class_mat An optional numeric matrix defining the color 
#' classes of each cell in the matrix. If omitted this is calculated
#' by the function.
#' @param grad_dir A string specifying how the color gradient 
#' should be applied. If "high" higher values are given darker
#' colors. If "low", lower values are given darker colors.
#' If "middle" values in the middle of the spectrum are 
#' given darker colors, and if "end" values at the ends
#' of the spectrum are given darker colors.
#' @param color_fun Either "linear" or "exponential" indicating
#' how the colors should transition from light to dark across
#' values.
#' @param exp_steepness If color_fun is "exponential," exp_steepness
#' indicates how quickly the colors should transition from light 
#' to dark.
#' @param global_color_scale Whether to impose a global minimum
#' and maximum to the colors, or to use the values themselves to 
#' determine the top and bottom of the color scale.
#' @param global_min If global_color_scale is TRUE, the minimum
#' value that should be assigned a color.
#' @param global_max If global_color_scale is TRUE, the maximum
#' value that should be assigned a color.
#' @param sig_digs The number of significant figures to use
#' from the input matrix. Helpful primarily when show_text is TRUE.
#' @param use_pheatmap_colors If TRUE, all other color parameters
#' are ignored, and colors like those used in the R package
#' pheatmap are used instead.
#' @param na_col The color to use for missing values.
#' @param gridlines Whether to plot gridlines on the matrix.
#'
#' @return None
#' 
#' @import pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics rasterImage
#' @keywords internal

image_with_text <- function(mat, xlab = "", ylab = "", main = NULL, main_shift = 0.12, 
	col_names = colnames(mat), row_names = rownames(mat), row_text_adj = 1, row_text_shift = 0, 
	row_text_rotation = 0, col_text_rotation = 90, col_text_adj = 1, 
	col_text_shift = 0, show_text = TRUE, cex = 0.5, col_text_cex = 1, 
	row_text_cex = 1, main_cex = 1, split_at_vals = FALSE, split_points = 0, 
	col_scale = "gray", color_spread = 50, light_dark = "f", class_mat = NULL, 
	grad_dir = c("high", "low", "middle", "ends"), color_fun = c("linear", "exponential"), 
	exp_steepness = 1, global_color_scale = FALSE, global_min = NULL, global_max = NULL, 
	sig_digs = 3, use_pheatmap_colors = FALSE, na_col = "lightgray", gridlines = FALSE){
		
		oldPar <- par(no.readonly = TRUE)
		on.exit(oldPar)
		
		#make sure Inf and -Inf are coded as NA
		mat[which(!is.finite(mat))] <- NA
		
		if(length(which(is.na(mat))) == length(mat)){
			return()
		}
		# if(length(light_dark) < length(col_scale)){light_dark <- rep(light_dark, length(col_scale))}
		
		get_default_col_fun <- grep("lin", color_fun)
		if(length(get_default_col_fun) > 0){
			color_fun = "linear"
		}
		
		end_fudge_factor = 10^-10

		if(is.null(class_mat)){
			class_mat <- matrix(1, dim(mat)[1], dim(mat)[2])
		}

		if(split_at_vals){
			for(p in 1:length(split_points)){
				class_mat[which(mat >= split_points[p])] <- class_mat[which(mat >= split_points[p])] + 1
			}
			# if(length(grad_dir) == 2){grad_dir <- "ends"}else{grad_dir <- "high"}		
		}else{
			split_points <- NULL	
		}
		class_mat[which(is.na(mat))] <- NA

		while(min(class_mat, na.rm = TRUE) > 1){class_mat <- class_mat - 1}
		
		classes <- sort(unique(as.vector(class_mat)))
		num_classes <- length(classes)
		if(num_classes == 1){
			class_mat <- matrix(1, nrow(mat), ncol(mat))
		}
	
			
		if(length(col_scale) == (length(split_points)+1)){
			class_cols <- col_scale
		}else{
			class_cols <- col_scale[classes]
		}
		if(length(col_scale) < num_classes){
			extra_cols <- num_classes - length(col_scale)
			col_scale <- c(col_scale, rep(col_scale, ceiling(extra_cols/length(col_scale))))
		}
			
		get_default <- grep("h", grad_dir)
		if(length(get_default) > 0){
			grad_dir <- "high"
		}
		
		# if(light_dark == "f"){max_col = 4}else{max_col = 8}
		max_col = 4
		dir_list <- vector(mode = "list", length = num_classes)
		names(dir_list) <- classes
		if(grad_dir == "high"){
			for(i in 1:length(dir_list)){
				dir_list[[i]] <- 1:max_col
			}
		}
		if(grad_dir == "low"){
			for(i in 1:length(dir_list)){
				dir_list[[i]] <- max_col:1
			}
		}
		if(grad_dir == "middle"){
			if(length(dir_list) != 2){stop("I can only color the middle if there are exactly two classes")}
			dir_list[[1]] <- 1:max_col
			dir_list[[2]] <- max_col:1
		}
			
		if(grad_dir == "ends"){
			if(length(dir_list) != 2){stop("I can only color the ends if there are exactly two classes")}
			dir_list[[1]] <- max_col:1
			dir_list[[2]] <- 1:max_col
		}

		#============================================================================
		#internal functions
		#============================================================================
		#This function takes in a matrix of values matched with colors, and 
		#a vector of values. It matched up the appropriate color for each
		#value in the vector
		bin_cols <- function(color_key, V){
			color_v <- rep(NA, length(V))
			for(i in 1:length(V)){
				diff_v  <- V[i] - as.numeric(color_key[,1])
				closest_val <- which(abs(diff_v) == min(abs(diff_v)))[1]
				color_v[i] <- color_key[closest_val,2]
			}
			return(color_v)
		}

		#This function generates the matrix of colors to use in 
		#the raster function
		fill_color_ramp <- function(mat, class_mat, global){

			#make color scales for each class or globally as defined
			color_scales <- vector(mode = "list", length = num_classes)
			names(color_scales) <- classes

			color_ramp <- matrix(NA, dim(mat)[1], dim(mat)[2])
			num_classes = length(unique(as.vector(class_mat)))
			
			for(cl in 1:length(classes)){
				if(global_color_scale){
					if(is.null(global_min)){min_cl = min(mat, na.rm = TRUE)}else{min_cl = global_min}
					if(is.null(global_max)){max_cl = max(mat, na.rm = TRUE)}else{max_cl = global_max}	
				}else{
					if(length(which(class_mat == classes[cl])) > 0 && !all(is.na(mat[which(class_mat == classes[cl])]))){
						min_cl <- min(mat[which(class_mat == classes[cl])], na.rm = TRUE)
						max_cl <- max(mat[which(class_mat == classes[cl])], na.rm = TRUE)
					}else{
						min_cl <- NA
					}
				}
					
				if(!is.na(min_cl)){
					if(color_fun == "linear"){
						color_levels <- seq(min_cl, max_cl, length=256)
					}else{
						color_levels <- exp_color_fun(min_cl, max_cl, steepness = exp_steepness, num_cols=256)	
					}
	
					#make the function to generate 
					col_vals <- get_color2(col_scale[cl], col_gap = color_spread)
										
					color_locale <- which(names(color_scales) == classes[cl])
					color_scales[[color_locale]] <- colorRampPalette(col_vals[dir_list[[color_locale]]])
					
					#find the entries in each class
					entry_locale <- which(class_mat == classes[cl])
					if(length(entry_locale) > 0){
						entry_vals <- mat[entry_locale]
						color_key <- cbind(color_levels, do.call(color_scales[[color_locale]], list(256)))
						entry_cols <- bin_cols(color_key, entry_vals)
						entry_order <- order(entry_locale)
						color_ramp[entry_locale[entry_order]] <- entry_cols
					}
				}
			}
			
			return(color_ramp)
		}
		
		#============================================================================


		if(use_pheatmap_colors){
			pal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
			if(global_color_scale){
				bks <- pheatmap_generate_breaks(c(global_min, global_max), length(pal), 
				center = F)
			}else{
				bks <- pheatmap_generate_breaks(mat, length(pal), center = F)
			}
			color_ramp <- pheatmap_scale_colours(mat, col=pal, breaks=bks, na_col = na_col)
		}else{
			color_ramp = fill_color_ramp(mat, class_mat, global_color_scale)
		}
		 	
		zmin <- min(mat, na.rm = TRUE)
		zmax <- max(mat, na.rm = TRUE)

		na_locale <- which(is.na(mat))
		if(length(na_locale) > 0){
			mat[na_locale] <- 0
		}
		
		all_ind <- which(!is.na(mat), arr.ind = TRUE)
		if(length(na_locale) > 0){
			mat[na_locale] <- NA
		}
		
		#translate the color_ramp matrix to rgb matrices so we can use grid.
		rgb_mat <- col2rgb(color_ramp)
		col <- matrix(rgb(rgb_mat["red",]/256,rgb_mat["green",]/256,rgb_mat["blue",]/256), dim(mat)[1], dim(mat)[2])

		max_dim <- max(c(dim(mat)[1], dim(mat)[2]))

		plot(c(1, dim(mat)[2]), c(1, dim(mat)[1]), type = "n", axes = FALSE, xlab = xlab, ylab = ylab, xlim = c(0.7, dim(mat)[2]+0.2), ylim = c(0.7, dim(mat)[1]+0.2), bg = "transparent")

		rasterImage(col, xleft = 0.5, ybottom = 0.5, xright = dim(mat)[2]+0.5, ytop = dim(mat)[1]+0.5, interpolate = FALSE, bg = "transparent")
		x_coord <- matrix(segment_region(0.5, dim(mat)[2]+0.5, dim(mat)[2], "center"), nrow = dim(mat)[1], ncol = dim(mat)[2], byrow = TRUE)
		y_coord <- matrix(segment_region(0.5, dim(mat)[1]+0.5, dim(mat)[1], "center"), nrow = dim(mat)[1], ncol = dim(mat)[2])
		
		if(show_text){
			text(x_coord, rev(y_coord), labels = signif(as.vector(mat), sig_digs), cex = cex)
		}
		
		if(!is.null(main)){
			par(xpd = TRUE)
			plot_range = max(y_coord) - min(y_coord)
			text(mean(x_coord[1,]), max(y_coord)+(plot_range*main_shift), labels = main, cex = main_cex)
			par(xpd = FALSE)
		}
			
		if(!is.null(col_names)){
			if(length(col_names) != dim(mat)[2]){
				stop("There is a different number of column names than columns.")
			}
			par(xpd = TRUE)
			# text(x_coord[1,], (min(y_coord)-(max(y_coord)) - col_text_shift), labels = col_names, srt = col_text_rotation, adj = col_text_adj)
			plot_range <- max(y_coord) - min(y_coord)
			text((x_coord[1,]), (min(y_coord) - (plot_range*0.1) - col_text_shift), labels = col_names, srt = col_text_rotation, adj = col_text_adj, cex = col_text_cex)
		}

		if(!is.null(row_names)){
			if(length(row_names) != dim(mat)[1]){
				stop("There is a different number of row names than rows.")
			}
			par(xpd = TRUE)
			plot_range <- (max(x_coord) - min(x_coord))
			text((min(x_coord) - (plot_range*0.1) - row_text_shift), y_coord[,1], labels = rev(row_names), adj = row_text_adj, cex = row_text_cex, srt = row_text_rotation)
		}


	}