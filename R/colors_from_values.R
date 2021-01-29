#' Retrieve colors based on numeric values
#' 
#' This function gets colors for numeric values
#' It can use one the function \code{\link{get_color}},
#' which accepts any of c("green", "purple", "orange", 
#' "blue", "brown", "gray"), or it can use pheatmap
#' colors.
#' The colors can be split into multiple classes based
#' on split_points. For example positive and negative 
#' numbers can each be plotted on their own scale, by 
#' setting split_at_vals to TRUE and setting split_points
#' to 0. split_points can be a vector defining multiple 
#' value classes for different colors. By default, the
#' different classes use the col_scale colors in order.
#' Specific colors for classes can be set using the 
#' col_scale argument.
#' 
#' @param vals A vector of numerical values for which to
#' retrieve colors.
#' @param split_at_vals A logical value indicating whether
#' the numbers should be broken into classes with different 
#' colors.
#' @param split_points If split_at_vals is TRUE, split_points
#' is used to define multiple classes of numbers. For example,
#' if split_points is 0, negative numbers will be assigned one
#' class of colors, and positive numbers will be assigned another.
#' @param col_scale One of c("green", "purple", "orange", "blue", 
#' "brown", "gray") to indicate the color scale to be used. Defaults
#' to gray.
#' @param light_dark One of "l", "d", or "f" indicating whether
#' the colors used should be light ("l"), dark ("d"), or from
#' across the full spectrum ("f").
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
#' @param use_pheatmap_colors If TRUE, all other color parameters
#' are ignored, and colors like those used in the R package
#' pheatmap are used instead.
#' @param na_col The color to use for missing values.
#'
#' @references Raivo Kolde (2019). pheatmap: Pretty Heatmaps. 
#' R package version 1.0.12. https://CRAN.R-project.org/package=pheatmap
#'
#' @return A vector of colors assigned to the number is vals is returned.
#' 
#' @import pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices col2rgb colorRampPalette rgb
#' @keywords internal

colors_from_values <- function(vals, split_at_vals = FALSE, split_points = 0, 
col_scale = c("green", "purple", "orange", "blue", "brown", "gray"), light_dark = "f", 
grad_dir = c("high", "low", "middle", "ends"), color_fun = c("linear", "exponential"), 
exp_steepness = 1, global_color_scale = FALSE, global_min = NULL, global_max = NULL, 
use_pheatmap_colors = FALSE, na_col = "lightgray"){

	 	class_mat = NULL
		#make sure Inf and -Inf are coded as NA
		vals[which(!is.finite(vals))] <- NA
		
		if(length(which(is.na(vals))) == length(vals)){
			return()
			}

		
		get_default_col_fun <- grep("lin", color_fun)
		if(length(get_default_col_fun) > 0){
			color_fun = "linear"
			}
		
		end_fudge_factor = 10^-10

		if(is.null(class_mat)){
			class_mat <- rep(1, length(vals))
			}

		if(split_at_vals){
			for(p in 1:length(split_points)){
				class_mat[which(vals >= split_points[p])] <- class_mat[which(vals >= split_points[p])] + 1
				}
			}else{
			split_points <- NULL	
			}
		class_mat[which(is.na(vals))] <- NA

		while(min(class_mat, na.rm = TRUE) > 1){class_mat <- class_mat - 1}
		
		classes <- sort(unique(as.vector(class_mat)))
		num_classes <- length(classes)
		if(num_classes == 1){
			class_mat <- rep(1, length(vals))
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
		fill_color_ramp <- function(vals, class_mat, global){

			#make color scales for each class or globally as defined
			color_scales <- vector(mode = "list", length = num_classes)
			names(color_scales) <- classes

			color_ramp <- rep(NA, length(vals))
			num_classes = length(unique(as.vector(class_mat)))
			
			for(cl in 1:length(classes)){
				if(global_color_scale){
					if(is.null(global_min)){min_cl = min(vals, na.rm = TRUE)}else{min_cl = global_min}
					if(is.null(global_max)){max_cl = max(vals, na.rm = TRUE)}else{max_cl = global_max}	
					}else{
					if(length(which(class_mat == classes[cl])) > 0 && !all(is.na(vals[which(class_mat == classes[cl])]))){
						min_cl <- min(vals[which(class_mat == classes[cl])], na.rm = TRUE)
						max_cl <- max(vals[which(class_mat == classes[cl])], na.rm = TRUE)
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
					col_vals <- get_color(col_scale[cl], light_dark)
					color_locale <- which(names(color_scales) == classes[cl])
					color_scales[[color_locale]] <- colorRampPalette(col_vals[dir_list[[color_locale]]])
					
					#find the entries in each class
					entry_locale <- which(class_mat == classes[cl])
					if(length(entry_locale) > 0){
						entry_vals <- vals[entry_locale]
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
				bks <- pheatmap_generate_breaks(vals, length(pal), center = F)
			}
			color_ramp <- pheatmap_scale_colours(vals, col=pal, breaks=bks, na_col = na_col)
		}else{
			color_ramp = fill_color_ramp(vals, class_mat, global_color_scale)
		}
		 	
		zmin <- min(vals, na.rm = TRUE); zmax <- max(vals, na.rm = TRUE)

		na_locale <- which(is.na(vals))
		if(length(na_locale) > 0){
			vals[na_locale] <- 0
			}
		
		all_ind <- which(!is.na(vals), arr.ind = TRUE)
		if(length(na_locale) > 0){
			vals[na_locale] <- NA
			}
		
		#translate the color_ramp matrix to rgb matrices so we can use grid.
		rgb_mat <- col2rgb(color_ramp)
		col <- rgb(rgb_mat["red",]/256,rgb_mat["green",]/256,rgb_mat["blue",]/256)

		return(col)
		

	}