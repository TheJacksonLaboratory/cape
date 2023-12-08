#' Plot trait pairs against each other
#' 
#' This function plots pairs of traits against each other
#' to visualize the correlations between traits.
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param pheno_which A vector of trait names to plot. The default is to plot all traits.
#' @param color_by A character string indicating a value to color the traits by, for example sex or treatment.
#' It must be one of the covariates. See \code{\link{pheno2covar}}.
#' @param group_labels A vector of names for the legend indicating the groups for the colored dots.
#' @param text_cex A numeric value indicating the size of the text
#' @param pheno_labels A vector of names for traits to appear in the plot in case the column names are not very pretty.
#' @param pt_cex A numeric value indicating the size of the points.
#' 
#' @importFrom graphics hist legend pairs rect
#' @importFrom stats cor
#' 
#' @export

plot_pheno_cor <- function(data_obj, pheno_which = NULL, color_by = NULL, group_labels = NULL,
  text_cex = 1, pheno_labels = NULL, pt_cex = 1){

	oldPar <- par(no.readonly = TRUE)
	on.exit(oldPar)
		
	all_pheno <- data_obj$pheno

	if(is.null(pheno_which)){
		pheno_names <- colnames(all_pheno)
	}else{
		if(is.numeric(pheno_which)[1]){
			pheno_names <- colnames(all_pheno)[pheno_which]
		}else{
			pheno_names <- pheno_which	
		}
	}
	
	sub_pheno <- all_pheno[,match(pheno_names, colnames(all_pheno))]
	
	if(is.null(pheno_labels)){
		pheno_labels <- pheno_names
	}
		
	pairs_matrix <- pair_matrix(pheno_names)

	all_cols <- brewer.pal(8, "Accent")

	if(!is.null(color_by)){
		group_mem <- NULL
		pheno_col <- which(colnames(data_obj$pheno) == color_by)
		if(length(pheno_col) > 0){
			group_mem <- data_obj$pheno[,pheno_col]
		}else{
			covar_info <- get_covar(data_obj)
			covar_col <- which(covar_info$covar_names %in% color_by)
			group_mem <- covar_info$covar_table[,covar_col]
		}			

		if(length(group_mem) == 0){
			stop(paste0("I couldn't find ", color_by, ". Please check the case and the spelling."))
		}
		
		cols <- rep(NA, dim(sub_pheno)[1])
		groups <- sort(unique(group_mem))
		num_groups <- length(groups)
		if(num_groups > 8){
			stop("There cannot be more than 8 groups")
		}
		for(i in 1:num_groups){
			cols[which(group_mem == groups[i])] <- all_cols[i]
		}
	}else{
		cols <- rep("black", dim(sub_pheno)[1])
		groups <- NULL
	}

	if(is.null(group_labels)){
		group_labels <- groups
	}

	plot_cor <- function(x,y){
		points(x,y, col = cols, cex = pt_cex, pch = 16)
		if(!is.null(color_by)){
			legend("topleft", legend = group_labels, fill = all_cols[1:length(group_labels)], cex = 0.7)
		}
	}
	
	write_cor <- function(x,y){
		if(length(which(!is.na(x+y))) > 0){
			if(!is.null(groups)){
				group_cor <- try(apply(matrix(groups, ncol = 1), 1, function(a) cor(x[which(group_mem == a)], y[which(group_mem == a)], use = "complete")), silent = TRUE)
			}else{
				group_cor = NULL
			}
			pheno_cor <- try(cor(x, y, use = "complete"), silent = TRUE)
		}

		cors <- c(group_cor, pheno_cor)
		group_labels <- c(group_labels, "R")
		x_shrinkage = 0.4
		y_shrinkage = 0.3
		x_range <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
		y_range <- max(y, na.rm = TRUE) - min(y, na.rm = TRUE)
		x_shift <- x_shrinkage*x_range
		y_shift <- y_shrinkage*y_range
		y_placement <- segment_region(min(y, na.rm = TRUE)+y_shift, max(y, na.rm = TRUE)-y_shift, max(c(length(cors),2)), alignment = "center")[1:length(cors)]
		x_placement <- segment_region(min(x, na.rm = TRUE)+x_shift, max(x, na.rm = TRUE)-x_shift,3, alignment = "center")
		text(rep(x_placement[1], length(cors)), y_placement, group_labels, cex = text_cex, adj = 1)
		text(rep(x_placement[2], length(cors)), y_placement, rep("=", length(cors)), cex = text_cex, adj = 0.5)
		text(rep(x_placement[3], length(cors)), y_placement, signif(cors, 3), cex = text_cex, adj = 0)	
	}
		
	panel_hist <- function(x, ...){
   		usr <- par("usr"); on.exit(par("usr" = usr))
    	par(usr = c(usr[1:2], 0, 1.5) )
	   	h <- hist(x, plot = FALSE)
    	breaks <- h$breaks; nB <- length(breaks)
	   	y <- h$counts; y <- y/max(y)
    	rect(breaks[-nB], 0, breaks[-1], y, ...)
	}

	pairs(sub_pheno, lower.panel = plot_cor, upper.panel = write_cor, diag.panel = panel_hist, labels = pheno_labels)
	# legend("topleft", legend = group_labels, fill = all_cols[1:length(group_labels)])

}
