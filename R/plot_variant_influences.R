#' Plot cape coefficients
#' 
#' This function plots the the cape coefficients between
#' pairs of markers as a heat map.
#' The interactions are shown in the main part of the heatmap
#' while the main effects are shown on the right hand side.
#' Directed interactions are read from the y axis to the x axis.
#' For example an interaction from marker1 to marker2 will be shown
#' in the row corresponding to marker1 and the column corresponding
#' to marker2. 
#' Similarly, if marker1 has a main effect on any traits, these
#' will be shown in the row for marker1 and the trait columns.
#' 
#' @param data_obj a \code{\link{Cape}} object
#' @param p_or_q A threshold indicating the maximum p value (or q value
#' if FDR was used) of significant interactions and main effects
#' @param min_std_effect An optional filter. The plot will exclude
#' all pairs with standardized effects below the number set here.
#' @param plot_all_vals If TRUE will plot all values regardless of 
#' significant
#' @param standardize Whether to plot effect sizes (FALSE) or standardized
#' effect sizes (TRUE)
#' @param color_scheme A character value of either "DO/CC" or other indicating the 
#' color scheme of main effects. If "DO/CC" allele effects can be plotted with the
#' DO/CC colors.
#' @param pos_col The color to use for positive main effects and interactions
#' must be one of "green", "purple", "red", "orange", "blue", "brown", "yellow", "gray"
#' see \code{\link{get_color}}
#' @param neg_col The color to use for negative main effects and interactions
#' takes the same values as pos_col.
#' @param not_tested_col The color to use for marker pairs not tested. Takes
#' the same values as pos_col and neg_col
#' @param show_marker_labels Whether to write the marker labels on the plot
#' @param show_chr Whether to show chromosome boundaries
#' @param label_chr Whether to label chromosomes if plotted
#' @param show_alleles If TRUE, the allele of each marker is indicated by color.
#' @param scale_effects One of "log10", "sqrt", "none." If some effects are
#' very large, scaling them can help show contrasts between smaller values.
#' The default is no scaling.
#' @param pheno_width Each marker and trait gets one column in the matrix. 
#' If there are many markers, this makes the effects on the traits difficult 
#' to see. pheno_width increases the number of columns given to the phenotypes.
#' For example, if pheno_width = 11, the phenotypes will be shown 11 times wider
#' than individual markers.
#' @param covar_width See pheno_width. This is the same effect for covariates.
#' @param covar_labels Labels for covariates if different from those stored in 
#' the data object.
#' @param phenotype_labels Labels for traits if different from those stored in 
#' the data object
#' @param show_not_tested Whether to color the marker pairs that were not
#' tested. If FALSE, they will not be colored in.
#' 
#' @return This function invisibly returns the variant influences matrix.
#' shown in the heat map.
#' 
#' @export
plot_variant_influences <- function(data_obj, p_or_q = 0.05, min_std_effect = 0, 
  plot_all_vals = FALSE, standardize = FALSE, 
  color_scheme = c("DO/CC", "other"), pos_col = "brown", neg_col = "blue", 
  not_tested_col = "lightgray", show_marker_labels = FALSE, show_chr = TRUE, 
  label_chr = TRUE, show_alleles = TRUE, scale_effects = c("log10", "sqrt", "none"), 
  pheno_width = NULL, covar_width = NULL, covar_labels = NULL, phenotype_labels = NULL, 
  show_not_tested = TRUE){
  
  if(!show_not_tested){
    not_tested_col = FALSE
  }
  
  geno_names <- data_obj$geno_names
  marker_names <- geno_names[[3]]
  
  if(length(grep("n", scale_effects)) > 0){
    scale_effects <- "none"
  }
  if(length(scale_effects) == 1){
    if(scale_effects != "log10" & scale_effects != "sqrt" & scale_effects != "none"){
      stop("scale_effects must be 'log10', 'sqrt' or 'none.'")
    }
  }
  
  var_influences <- data_obj$var_to_var_p_val
  
  pheno_inf <- data_obj$max_var_to_pheno_influence
  if(is.null(phenotype_labels)){
    pheno_names <- names(data_obj$max_var_to_pheno_influence)
  }else{
    pheno_names <- phenotype_labels
    if(length(pheno_names) != length(names(data_obj$max_var_to_pheno_influence))){
      stop("I am detecting the wrong number of phenotype labels for the phenotypes present.")
    }
  }
  num_pheno <- length(pheno_names)
  
  if(not_tested_col == TRUE){
    not_tested_col = "lightgray"
  }
  
  if(is.null(var_influences)){
    stop("calc_p() must be run to calculate variant-to-variant influences.")
  }
  
  if(is.null(pheno_inf)){
    stop("direct_influence() must be run to calculate variant-to-trait influences.")
  }
  
  #This function expands the given rows or columns of 
  #a matrix by a given amount
  expand_matrix <- function(mat, row_col_num, row_or_col, expansion_factor){
    if(expansion_factor == 1){
      return(mat)
      }
    if(row_or_col == "row"){
      row_labels <- 1:nrow(mat)
      for(i in 1:length(row_col_num)){
        mat_before <- mat[which(row_labels < row_col_num[i]),,drop=FALSE]
        mat_after <- mat[which(row_labels > row_col_num[i]),,drop=FALSE]
        row_to_expand <- mat[which(row_labels == row_col_num[i]),,drop=FALSE]
        mat_to_add <- matrix(row_to_expand, nrow = expansion_factor, ncol = ncol(mat), byrow = TRUE)
        label_locale <- which(row_labels == row_col_num[i])
        label <- rownames(mat)[label_locale]
        rowname_v <- rep("", expansion_factor)
        rowname_v[round(expansion_factor/2)] <- label
        rownames(mat_to_add) <- rowname_v
        row_labels <- c(row_labels[which(row_labels < row_col_num[i])], rep(row_col_num[i], expansion_factor), row_labels[which(row_labels > row_col_num[i])])
        mat <- rbind(mat_before, mat_to_add, mat_after)
      }
      return(mat)
    }
    
    
    if(row_or_col == "col"){
      col_labels <- 1:ncol(mat)
      for(i in 1:length(row_col_num)){
        mat_before <- mat[,which(col_labels < row_col_num[i]),drop=FALSE]
        mat_after <- mat[,which(col_labels > row_col_num[i]),drop=FALSE]
        col_to_expand <- mat[,which(col_labels == row_col_num[i]),drop=FALSE]
        mat_to_add <- matrix(col_to_expand, ncol = expansion_factor, nrow = nrow(mat), byrow = FALSE)
        
        label_locale <- which(col_labels == row_col_num[i])
        label <- colnames(mat)[label_locale]
        colname_v <- rep("", expansion_factor)
        colname_v[round(expansion_factor/2)] <- label
        colnames(mat_to_add) <- colname_v
        
        col_labels <- c(col_labels[which(col_labels < row_col_num[i])], rep(row_col_num[i], expansion_factor), col_labels[which(col_labels > row_col_num[i])])
        mat <- cbind(mat_before, mat_to_add, mat_after)
      }
      return(mat)
    }
  }
  
  
  unique_markers <- unique(c(as.vector(var_influences[,"Source"]), as.vector(var_influences[,"Target"]), rownames(pheno_inf[[1]])))
  split_markers <- strsplit(unique_markers, "_")
  just_markers <- sapply(split_markers, function(x) x[1])
  just_alleles <- sapply(split_markers, function(x) x[2])
  unique_marker_locale <- match(just_markers, marker_names)
  marker_order <- order(unique_marker_locale)		
  sorted_markers <- unique_markers[marker_order]
  
  if(show_alleles){
    allele_colors <- get_allele_colors(color_scheme, just_alleles)
    allele_cols <- allele_colors[match(just_alleles[marker_order], allele_colors[,2]),3]
  }else{
    allele_cols <- NULL
  }
  
   
  #update the markers based on the covariate width
  covar_info <- get_covar(data_obj)
  covar_names <- covar_info$covar_names
  if(is.null(covar_width)){
    covar_width <- round((length(sorted_markers)+length(covar_names))/20)
    if(covar_width < 1){
      covar_width = 1
    }
  }
  
  
  if(length(covar_names) > 0){
    covar_markers <- covar_names
    covar_locale <- match(covar_markers, sorted_markers)
    sorted_markers[sort(covar_locale)] <- covar_names #make sure the names are in the right order
    new_covar_markers <- sort(rep(covar_markers, covar_width))
    just_markers <- sorted_markers[-covar_locale]
    expanded_markers <- c(just_markers, new_covar_markers)
    #extend allele cols with gray for covariates
    allele_cols <- c(allele_cols, rep("lightgray", length(new_covar_markers)))
  }else{
    expanded_markers <- sorted_markers
    covar_locale <- NULL	
  }
  
  
  
  #get coordinates of the chromosome boundaries
  if(show_chr){
    num_covar <- length(covar_names)
    orig_chromosomes <- get_marker_chr(data_obj, markers =  sapply(strsplit(sorted_markers, "_"), function(x) x[[1]]))
    covar_locale <- which(orig_chromosomes == 0)
    chromosomes <- orig_chromosomes[which(orig_chromosomes != 0)]
    chr_labels <- rep(1, length(unique(chromosomes)))
    if(num_covar > 0){
      chr_labels <- c(chr_labels, rep(0, num_covar))
      for(i in 1:length(covar_names)){
        chromosomes <- c(chromosomes, rep(covar_names[i], covar_width))
      }
    }
    
    u_chr <- unique(chromosomes[which(!is.na(chromosomes))])
    chr_boundaries <- apply(matrix(u_chr, ncol = 1), 1, function(x) max(which(chromosomes == x))) + 0.5
    chr_boundaries <- c(0, chr_boundaries)
    if(label_chr){
      #use only the first two characters of each chromosome
      chr_names <- unique(chromosomes)
      if(!is.null(covar_labels)){
        chr_names[(length(chr_names)-num_covar+1):length(chr_names)] <- covar_labels
      }
      
      # chr_names <- unlist(lapply(strsplit(unique(chromosomes), ""), function(x) paste(x[1:2], collapse = "")))
      # chr_names <- unlist(strsplit(chr_names, "NA"))
    }else{
      chr_names <- NULL
    }
  }else{
    chr_boundaries <- NULL
    chr_names <- NULL
  }
  
  #make the variant influence matrix from significant interactions
  var_sig_col <- which(colnames(var_influences) == "p_adjusted")
  
  if(standardize){
    edge_weights <- as.numeric(var_influences[,5])
  }else{
    edge_weights <- as.numeric(var_influences[,3])	
  }
  pvals <- as.numeric(var_influences[,var_sig_col])
  
  var_influence_net <- graph_from_edgelist(var_influences[,1:2], directed = TRUE)
  var_pval_net <- var_influence_net
  E(var_influence_net)$weight <- edge_weights
  E(var_pval_net)$weight <- pvals
  var_influence_mat <- as.matrix(as_adjacency_matrix(var_influence_net, attr = "weight"))
  var_pval_mat <- as.matrix(as_adjacency_matrix(var_pval_net, attr = "weight"))
  
  #put in marker order
  mat_order <- match(sorted_markers, rownames(var_influence_mat))
  var_influence_mat <- var_influence_mat[mat_order,mat_order]
  var_pval_mat <- var_pval_mat[mat_order,mat_order]

  #turn the not-tested elements to NA
  var_influence_mat[which(var_influence_mat == 0)] <- NA
  var_pval_mat[which(var_pval_mat == 0)] <- NA
  
  pheno_influence_mat <- matrix(NA, nrow = length(unique_markers), ncol = num_pheno)
  pheno_pval_mat <- matrix(NA, nrow = length(unique_markers), ncol = num_pheno)
  colnames(pheno_influence_mat) <- colnames(pheno_pval_mat) <- pheno_names
  rownames(pheno_influence_mat) <- rownames(pheno_pval_mat) <- sorted_markers
  
  
  #expand the covariate rows and columns to the specified width
  if(length(covar_names) > 0){
    covar_locale <- which(sorted_markers %in% covar_names)
    
    new_var_inf <- expand_matrix(mat = var_influence_mat, 
      row_col_num = covar_locale, row_or_col = "row", 
      expansion_factor = covar_width)
    new_var_inf <- expand_matrix(new_var_inf, covar_locale, "col", covar_width)
    
    new_var_pval <- expand_matrix(var_pval_mat, covar_locale, "row", covar_width)
    new_var_pval <- expand_matrix(new_var_pval, covar_locale, "col", covar_width)		
    
    var_influence_mat <- new_var_inf
    var_pval_mat <- new_var_pval
  }
  
  
  
  #fill the variant-to-phenotype matrix with test statistics 
  #(still with sources in rows and targets in columns)
  #use phenotypes or eigentraits based on user input
  pheno_sig_col <- which(colnames(pheno_inf[[1]]) == "p_adjusted")
  for(i in 1:length(unique_markers)){
    for(j in 1:length(pheno_names)){
      marker_locale <- which(pheno_inf[[j]][,1] == unique_markers[i])
      if(length(marker_locale) > 0){
        if(standardize){	
          pheno_influence_mat[as.character(unique_markers[i]), pheno_names[j]] <- pheno_inf[[j]][marker_locale, "t_stat"]
        }else{
          pheno_influence_mat[as.character(unique_markers[i]), pheno_names[j]] <- pheno_inf[[j]][marker_locale, "coef"]
        }
        pheno_pval_mat[as.character(unique_markers[i]), pheno_names[j]] <- pheno_inf[[j]][marker_locale, pheno_sig_col]
      }else{
        pheno_influence_mat[as.character(unique_markers[i]), pheno_names[j]] <- NA
        pheno_pval_mat[unique_markers[i], pheno_names[j]] <- NA
      }
    }
  }
  
  if(is.null(pheno_width)){
    pheno_width <- round((length(sorted_markers)+length(covar_names))/15)
    if(pheno_width < 1){
      pheno_width <- 1
      }
  }
  #expand the phenotype influence matrix to give it more visual weight in the plot
  expanded_pheno_mat <- expand_matrix(mat = pheno_influence_mat, 1:ncol(pheno_influence_mat), "col", pheno_width)
  expanded_pheno_pval_mat <- expand_matrix(pheno_pval_mat, 1:ncol(pheno_pval_mat), "col", pheno_width)
  
  #also expand the regions where the covariates are if there are covariates
  if(length(covar_names) > 0){
    expanded_pheno_mat <- expand_matrix(expanded_pheno_mat, covar_locale, "row", covar_width)
    expanded_pheno_pval_mat <- expand_matrix(expanded_pheno_pval_mat, covar_locale, "row", covar_width)
  }
  
  
  full_inf_mat <- cbind(var_influence_mat, expanded_pheno_mat)
  full_pval_mat <- cbind(var_pval_mat, expanded_pheno_pval_mat) 
  
  full_inf_mat_num <- apply(full_inf_mat, 2, as.numeric)
  rownames(full_inf_mat_num) <- rownames(full_inf_mat)
  colnames(full_inf_mat_num) <- colnames(full_inf_mat)
  
  full_pval_mat_num <- apply(full_pval_mat, 2, as.numeric)
  dimnames(full_pval_mat_num) <- dimnames(full_pval_mat)

  
  #get the coordinates for all pairs not tested
  # not_tested_locale <- which(is.na(rotate_mat(full_inf_mat)), arr.ind = TRUE)
  not_tested_locale <- which(is.na(rotate_mat(full_inf_mat_num)), arr.ind = TRUE)
  
  if(not_tested_col == FALSE || is.na(not_tested_col)){
    not_tested_locale <- NULL
  }
  
  #take out any values that aren't significant according
  #to the user cutoff, and do not have a high enough
  #effect size
  #if we are not plotting all value
  #use an extra color matrix to highlight significant interactions
  #if we are plotting all markers
  
  extra_col_mat <- NULL
  if(!plot_all_vals){
    full_inf_mat_num[which(full_pval_mat_num > p_or_q)] <- NA
    full_inf_mat_num[which(abs(full_inf_mat_num) < min_std_effect)] <- NA
  }else{
    extra_col_mat <- matrix(NA, nrow = nrow(full_pval_mat_num), ncol = ncol(full_pval_mat_num))
    
    min_effect <- which(abs(full_inf_mat_num) > min_std_effect)
    neg_effect <- intersect(min_effect, which(full_inf_mat_num < 0))
    pos_effect <- intersect(min_effect, which(full_inf_mat_num > 0))
    
    sig_neg <- intersect(neg_effect, which(full_pval_mat_num < p_or_q))
    sig_pos <- intersect(pos_effect, which(full_pval_mat_num < p_or_q))
    
    
    extra_col_mat[sig_neg] <- get_color(neg_col, light_dark = "d")[3]
    extra_col_mat[sig_pos] <- get_color(pos_col, light_dark = "d")[3]
  }
  main <- "Variant Influences"
  
  if(scale_effects == "log10"){
    neg_locale <- which(full_inf_mat_num < 0)
    scaled_effects <- log10(abs(full_inf_mat_num))
    scaled_effects[neg_locale] <- scaled_effects[neg_locale]*-1	
    full_inf_mat_num <- scaled_effects
    main <- "log10 Variant Influences"
  }
  if(scale_effects == "sqrt"){
    neg_locale <- which(full_inf_mat_num < 0)
    scaled_effects <- sqrt(abs(full_inf_mat_num))
    scaled_effects[neg_locale] <- scaled_effects[neg_locale]*-1	
    full_inf_mat_num <- scaled_effects
    main <- "Square Root of Variant Influences"
  }
  
  
  if(length(which(na.omit(as.vector(full_inf_mat_num)) != 0)) == 0){
    plot.new()
    plot.window(xlim = c(0,1), ylim = c(0,1))
    text(0.5, 0.5, "No Significant Interactions")
  }else{
    
    min_val <- min(full_inf_mat_num, na.rm = TRUE)
    max_val <- max(full_inf_mat_num, na.rm = TRUE)
    if(min_val == max_val){
      if(min_val < 0){max_val <- 0}
      if(min_val > 0){min_val <- 0}
    }
    my_image_plot(x = full_inf_mat_num, min_x = min_val, max_x = max_val, 
      main = main, xlab = "Target", ylab = "Source", 
      mark_coords = not_tested_locale, mark_col = not_tested_col, 
      show_labels = show_marker_labels, chromosome_coordinates = chr_boundaries, 
      chr_names = chr_names, chr_labels = chr_labels, 
      show_pheno_labels = TRUE, extra_col_mat = extra_col_mat, 
      allele_cols = allele_cols)
    
    #add phenotype names
    if(!is.null(not_tested_locale)){
      legend("topright", legend = "not testable", col = not_tested_col, pch = 16)
    }
  }
  
  invisible(full_inf_mat_num)
  
  
}