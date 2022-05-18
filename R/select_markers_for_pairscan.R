#' Select markers for the pairwise scan.
#'
#' This function selects markers for the pairwise scan.
#' Beause Cape is computationally intensive, pairscans 
#' should not be run on large numbers of markers. 
#' As a rule of thumb, 1500 markers in a population of 
#' 500 individuals takes about 24 hours to run without the
#' kinship correction. The kinship correction increases the
#' time of the analysis, and users may wish to reduce the number
#' of markers scanned even further to accommodate the extra
#' computational burden of the kinship correction.
#'
#' This function can select markers either from a pre-defined list
#' input as the argument \code{specific_markers}, or can select
#' markers based on their main effect size.
#' 
#' To select markers based on main effect size, this function 
#' first identifies effect score peaks using an automated
#' peak detection algorithm. It finds the peaks rising 
#' above a starting threshold and samples markers within each
#' peak based on the user-defined sampling density \code{peak_density}.
#' Setting \code{peak_density} to 0.5 will result in 50\% of the markers
#' in a given peak being sampled uniformly at random. Sampling
#' reduces the redundancy among linked markers tested in the pairscan.
#' If LD is relatively low in the population, this density can be
#' increased to 1 to include all markers under a peak. If LD is high,
#' the density can be decreased to reduce redundancy further. 
#' 
#' The algorithm compares the number of markers sampled to the target
#' defined by the user in the argument \code{num_alleles}. If fewer
#' than the target have been selected, the threshold is lowered, and 
#' the process is repeated until the target number of alleles have
#' been selected (plus or minus the number set in \code{tolerance}).
#' 
#' If the number of target alleles exceeds the number of markers
#' genotyped, all alleles will be selected automatically. 
#' 
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param singlescan_obj a singlescan object from \code{\link{singlescan}}.
#' @param geno_obj a genotype object
#' @param specific_markers A vector of marker names specifying which
#' markers should be selected. If NULL, the function uses main effect
#' size to select markers.
#' @param num_alleles The target number of markers to select if using
#' main effect size
#' @param peak_density The fraction of markers to select under each 
#' peak exceeding the current threshold. Should be set higher for populations
#' with low LD. And should be set lower for populations with high LD. Defaults 
#' to 0.5, corresponding to 50\% of markers selected under each peak.
#' @param window_size The number of markers to use in a smoothing window when
#' calculating main effect peaks. If NULL, the window size is selected automatically
#' based on the number of markers with consecutive rises and falls of main effect
#' size.
#' @param tolerance The allowable deviation from the target marker number in
#' number of markers. For example, If you ask the function to select 100 markers,
#' an set the tolerance to 5, the algorithm will stop when it has selected between
#' 95 and 105 markers.
#' @param plot_peaks Whether to plot the singlescan peaks identified by \code{\link{bin_curve}}.
#' This can be helpful in determining whether the window_size and peak_density parameters 
#' are optimal for the population.
#' @param verbose Whether progress should be printed to the screen
#' @param pdf_filename If plot_peaks is TRUE, this argument specifies the filename
#' to which the peaks are plotted.
#' 
#' @seealso \code{\link{bin_curve}}, \code{\link{singlescan}}
#' 
#' @return Returns the \code{\link{Cape}} object with a new matrix called
#' \code{geno_for_pairscan} containing the genotypes of the selected markers
#' for each individual.
#' 
#' @importFrom abind adrop
#' @importFrom grDevices dev.off pdf 
#' @importFrom graphics layout
#' @importFrom stats runif
#' 
#' @examples 
#' \dontrun{
#' #select 100 alleles to run through pairscan
#' data_obj <- select_markers_for_pairscan(data_obj, singlescan_obj, geno_obj, num_alleles = 100)
#' }
#' 
#' @export
select_markers_for_pairscan <- function(data_obj, singlescan_obj, geno_obj, 
  specific_markers = NULL, num_alleles = 50, peak_density = 0.5, window_size = NULL, 
  tolerance = 5, plot_peaks = FALSE, verbose = FALSE, pdf_filename = "Peak.Plots.pdf"){
  
  # These two lines need to be commented out on the VM, as of now, the pdf device is not supported
  oldPar <- par(no.readonly = TRUE)
  on.exit(oldPar)
  
  # If plot_pdf is FALSE we change the extension to .jpg
  if (endsWith(pdf_filename, '.pdf')) {
    if (!data_obj$plot_pdf) {
      pdf_filename <- gsub(".pdf",  ".jpg", pdf_filename)
    }
  }
  chr <- unique(data_obj$chromosome)
  
  geno <- get_geno(data_obj, geno_obj)
  alleles <- dimnames(geno)[[2]]
  n_alleles <- length(alleles)
  
  class_singlescan <- class(singlescan_obj)
  if(class_singlescan == "list"){
    ref_allele <- singlescan_obj$ref_allele
    data_obj$ref_allele <- ref_allele
  }else{
    ref_allele <- data_obj$ref_allele
  }
  
  if(is.null(ref_allele)){
    allele_text <- paste(alleles, collapse = ", ")
    ref_allele <- readline(prompt = paste("Which allele do you want to use as the reference?\n", allele_text, "\n"))
    data_obj$ref_allele <- ref_allele
  }
  
  # if both num_alleles is defined and the number is > the number of markers in the geno 
  # data, then ALL the markers are selected. However, we only want to allow this if
  # specific_markers is undefined. 
  
  if(is.null(specific_markers)) {
    #if we are asking for more markers than there are in the dataset
    #just take all of them.
    if (num_alleles >= dim(geno)[3]*(dim(geno)[2]-1))  {
      num_alleles = dim(geno)[3]
      alt_alleles <- setdiff(dimnames(geno)[[2]], ref_allele)
      specific_markers <- paste(dimnames(geno)[[3]], alt_alleles, sep = "_")
    }
  }
  
  if(!is.null(specific_markers)) {
    if(n_alleles == 2){
      ref_allele_locale <- which(data_obj$geno_names[[2]] == ref_allele)
      other_allele <- setdiff(1:2, ref_allele_locale)
      split_markers <- strsplit(as.character(specific_markers), "_")
      just_markers <- sapply(split_markers, function(x) x[1])
      just_marker_locale <- match(just_markers, dimnames(geno)[[3]])
      just_marker_locale <- just_marker_locale[which(!is.na(just_marker_locale))]
      geno_for_pairscan <- geno[,other_allele,just_marker_locale]
      colnames(geno_for_pairscan) <- paste(colnames(geno_for_pairscan), data_obj$geno_names[[2]][other_allele], sep = "_")
    }else{
      split_markers <- strsplit(specific_markers, "_")
      just_markers <- sapply(split_markers, function(x) x[1])
      just_alleles <- sapply(split_markers, function(x) x[2])
      geno_for_pairscan <- matrix(nrow = nrow(data_obj$pheno), ncol = length(just_markers))
      colnames(geno_for_pairscan) <- specific_markers
      for(i in 1:length(just_markers)){
        geno_for_pairscan[,i] <- geno[,just_alleles[i], just_markers[i]]
      }
    }
    
    if(verbose){cat("Removing markers that are not linearly independent...\n")}
    data_obj$geno_for_pairscan <- geno_for_pairscan
    geno_ind <- get_linearly_independent(data_obj)
    if(verbose){
      cat(length(geno_ind[[2]]), "allele(s) rejected.\n")
      cat("Final alleles selected:", "\t", ncol(geno_ind$independent_markers), "\n")
    }
    #we still need to specify values for selecting markers for the pairscan null distribution
    #so just use default value
    # cat("Generating the null distribution for the pairscan requires values for peak_density, window_size, and tolerance.\nSetting default values...")
    # data_obj$peak_density <- peak_density
    # data_obj$window_size <- window_size
    # data_obj$tolerance <- tolerance
    return(data_obj)
  }
  

  if(class_singlescan == "list"){ 
    results <- abs(singlescan_obj$singlescan_t_stats) #an actual singlescan object
  }else{
    results <- abs(singlescan_obj) #a singlescan matrix for calculating pairscan null distribution
  }
  
  
  filtered_results <- results
  
  covar_info <- get_covar(data_obj)
  results_no_covar <- results[which(!rownames(results) %in% covar_info$covar_names),,,drop=FALSE]
  result_chr <- get_marker_chr(data_obj, markers = rownames(results_no_covar), character_names = TRUE)  
    
  #===============================================================
  #internal functions
  #===============================================================
  
  #how may peaks are above a given cutoff?
  num_peaks <- function(allele_curves, bins, cutoff){
    filtered_bins <- bins
    filtered_bins[which(abs(allele_curves) < cutoff)] <- NA
    num_peaks <- apply(filtered_bins, 2, function(x) length(unique(x))-1)
    return(num_peaks)
  }
  
  #how many alleles are above a given t stat cutoff
  num_markers <- function(allele_curves, cutoff){
    filtered_curves <- allele_curves
    filtered_curves[which(abs(allele_curves) < cutoff)] <- NA
    num_markers <- apply(filtered_curves, 2, function(x) length(which(!is.na(x))))
    return(num_markers)
  }
  
  #how many markers are in each peak at a given cutoff?
  markers_per_peak <- function(allele_curves, bins, cutoff){
    #make a results mat with enough columns for each allele
    #and enough rows for each bin. Each cell will count the
    #number of markers in the bin at the designated cutoff
    result_mat <- matrix(0, ncol = dim(allele_curves)[[3]], nrow = max(bins, na.rm = TRUE))
    rownames(result_mat) <- 1:nrow(result_mat)
    colnames(result_mat) <- dimnames(allele_curves)[[3]]
    
    #delete all effects that are less than the cutoff
    filtered_bins <- bins
    filtered_bins[which(abs(allele_curves) < cutoff)] <- NA
    #image(filtered_bins)
    #count the number of markers in each bin for each allele
    for(i in 1:ncol(filtered_bins)){
      counts <- table(filtered_bins[,i])
      result_mat[names(counts),i] <- counts
    }
    return(result_mat)
  }

  #find how many markers in each peak if we sample at the 
  #specified density, with at least one marker per peak.
  num_sampled_markers  <- function(num_per_peak, peak_density){
    final_count <- num_per_peak
    above_thresh <- which(num_per_peak > 0) #find peaks with effect sizes above threshold
    final_count[above_thresh] <- 1 + floor(num_per_peak[above_thresh]*peak_density)
    return(final_count)
  }
  
  sample_peaks <- function(pheno_results, num_per_peak, bins){
    #sample markers based on peaks across all alleles
    sampled_markers <- vector(mode = "list", length = ncol(pheno_results))
    names(sampled_markers) <- colnames(pheno_results)
    for(i in 1:ncol(pheno_results)){
      #figure out which peaks we will sample from
      allele_markers <- NULL
      peaks_which <- which(num_per_peak[,i] > 0)
      if(length(peaks_which) > 0){
        for(j in 1:length(peaks_which)){
          #in each peak, pick the max, and sample the rest
          marker_locale <- which(bins[,i] == peaks_which[j])
          allele_markers <- c(allele_markers, marker_locale[which.max(abs(pheno_results[marker_locale,i]))])
          num_to_sample <- num_per_peak[peaks_which[j],i] - 1 #take off the maximum marker
          if(num_to_sample > 0){ #if there are still markers to get after grabbing the max
            unif_markers <- round(runif(num_to_sample, min = min(marker_locale), max = max(marker_locale)))
            allele_markers <- c(allele_markers, unif_markers)
          }#end case for sampling peak uniformly
        }#end looping through peaks for one allele
        sampled_markers[[i]] <- rownames(pheno_results)[sort(allele_markers)]
      }#end looping through alleles
    }
    return(sampled_markers)
  }
  
  #===============================================================
  
  #===============================================================
  #group markers into bins based on their effect size profiles
  #===============================================================	
  num_pheno <- dim(results)[[2]]
  allele_bins <- vector(mode = "list", length = num_pheno)
  for(ph in 1:num_pheno){
    if(verbose){cat("\nBinning markers for", colnames(filtered_results)[ph], "\n")}
    pheno_results <- results_no_covar[,ph,,drop=FALSE]
    
    if(plot_peaks){
      if (data_obj$plot_pdf) {
        cat(paste("Plotting results_no_covar: ", pdf_filename))
        pdf(pdf_filename, width = nrow(results_no_covar)*0.5, height = 15) 
      }
      cat(paste("Plotting results_no_covar: ", pdf_filename))
      jpeg(pdf_filename, res = 400, width = nrow(results_no_covar)*0.5, height = 15, units = "in")
      
      layout_mat <- get_layout_mat(ncol(pheno_results), "upright")
      # quartz(width = 15, height = 15)
      layout(layout_mat)
      par(mar = c(2,2,2,2))
    }
    
    #bin markers by chromosome
    for(ch in 1:length(chr)){
      if(verbose){
        report_progress(ch, length(chr))
      }

      chr_locale <- which(result_chr == chr[ch])
      chr_results <- pheno_results[chr_locale,,,drop=FALSE]
      
      #bin each allele curve by chromosome, so we don't have bins overlapping 
      #chromosome breaks. Make sure the first bin on one chromosome is 1+ the max
      #bin on the last chromosome
      chr_bins <- apply(chr_results, 3, function(x) bin_curve(x, plot_peaks = plot_peaks, window_size = window_size)$bins)
      if(!is.null(allele_bins[[ph]])){ 
        max_bin <- apply(allele_bins[[ph]], 2, max)
      }else{
        max_bin <- rep(0, ncol(chr_bins))#if we don't have any bins yet, initialize at 0
      }
      total_bins <- Reduce("cbind", lapply(1:ncol(chr_bins), function(x) chr_bins[,x]+max_bin[x]))
      total_bins <- matrix(total_bins, ncol = ncol(chr_bins))
      allele_bins[[ph]] <- rbind(allele_bins[[ph]], total_bins)
    }
  }
  if(verbose){cat("\n")}
  
  #===============================================================
  
  
  #===============================================================			
  #Determine the number of alleles to sample in each peak
  #lower the effect size cutoff until we have the approximate 
  #number of alleles requested
  #===============================================================
  if(verbose){cat("\nFinding effect size threshold...\n")}
  sorted_results <- sort(abs(as.vector(results_no_covar)))
  #start with a cutoff that might be near the number of 
  #alleles desired
  guess_point <- num_alleles
  min_effect_size = min(tail(sorted_results, guess_point))

  total_alleles <- 0
  alleles_checked <- NULL
  repeats <- 0 #checks for no change in allele number
  flips <- 0 #checks for bouncing around up and down around
  #around the target number
  step_size <- 0.1 #starting effect size step
  last_round <- "under"
  #while we are far away from the target number, and we have not
  #hit any of our checks.
  while((total_alleles < (num_alleles - tolerance) || total_alleles > (num_alleles + tolerance)) && repeats <= 100 && flips < 2){
    
    if(total_alleles < num_alleles){
      this_round  <- "under"
    }else{
      this_round <- "over"
    }

    if(this_round != last_round){
      flips <- flips + 1
    }

    #if we've flipped across the threshold, decrease the step size by 50%
    if(flips > 0){
      step_size <- step_size/2
    }

    ph_alleles <- vector(mode = "list", length = num_pheno)
    
    for(ph in 1:num_pheno){
      pheno_results <- results_no_covar[,ph,,drop=FALSE]
      ph_alleles[[ph]] <- markers_per_peak(allele_curves = pheno_results, 
      bins = allele_bins[[ph]], cutoff = min_effect_size)
    }

    num_sampled_alleles <- lapply(ph_alleles, function(x) num_sampled_markers(x, peak_density))
    total_alleles <- sum(unlist(num_sampled_alleles))

    if(verbose){
      cat(signif(min_effect_size,3), "\t\t", total_alleles, "\n")
    }
    
    #update and check number of repeats
    if(!is.null(alleles_checked) && tail(alleles_checked, 1) == total_alleles){
      repeats = repeats + 1
    }else{
      repeats <- 0 #otherwise, reset repeats
    }

    if(repeats > 1){
      #if we're starting to level off in the number of alleles we're
      #finding, increase the step.size by 50%
      step_size <- step_size *1.5
    }

    #update effect size threshold
    if(total_alleles < num_alleles){
      min_effect_size = min_effect_size - step_size
    }else{
      min_effect_size = min_effect_size + step_size
    }
    
    alleles_checked <- c(alleles_checked, total_alleles)
    
    last_round <- this_round
  }
  
  #sample markers from each peak as defined by num_sampled_alleles
  sampled_markers <- vector(mode = "list", length = num_pheno)
  for(ph in 1:num_pheno){	
    #phenotype results across all alleles
    #should be a matrix in two dimensions with 
    #individuals in rows and alleles in columns
    pheno_results <- results_no_covar[,ph,,drop=FALSE]
    pheno_results <- adrop(pheno_results, drop = 2)
    sampled_markers[[ph]] <- sample_peaks(pheno_results = pheno_results, 
    num_per_peak = num_sampled_alleles[[ph]], bins = allele_bins[[ph]])
  }
  
  #build a genotype matrix just with the sampled alleles
  num_parents <- ncol(pheno_results)
  markers_by_parent <- vector(mode = "list", length = num_parents)
  names(markers_by_parent) <- colnames(pheno_results)
  for(p in 1:num_parents){
    markers_by_parent[[p]] <- unique(unlist(lapply(sampled_markers, function(x) x[[colnames(pheno_results)[p]]])))
  }
  
  total_markers <- length(unlist(markers_by_parent))
  if(verbose){cat("total unique alleles:", "\t", total_markers, "\n")}
    
  geno_for_pairscan <- matrix(NA, nrow = nrow(data_obj$pheno), ncol = total_markers)
  colnames(geno_for_pairscan) <- 1:ncol(geno_for_pairscan)
  start_allele <- 1
  for(p in 1:length(markers_by_parent)){
    all_allele_markers <- markers_by_parent[[p]]
    if(length(all_allele_markers) > 0){
      geno_section <- geno[,colnames(pheno_results)[p],all_allele_markers,drop=FALSE]
      geno_for_pairscan[,start_allele:(start_allele+length(all_allele_markers)-1)] <- geno_section
      colnames(geno_for_pairscan)[start_allele:(start_allele+length(all_allele_markers)-1)] <- paste(all_allele_markers, colnames(pheno_results)[p], sep = "_")
      start_allele <- start_allele+length(all_allele_markers)
    }
  }
  
  #now order the matrix first by marker then by allele
  marker_split <- strsplit(colnames(geno_for_pairscan), "_")
  just_markers <- unlist(lapply(marker_split, function(x) x[1]))
  just_allele <- unlist(lapply(marker_split, function(x) x[2]))
  marker_table <- cbind(get_marker_num(data_obj, just_markers), just_allele)
  sorted_table <- sort_by_then_by(marker_table, col_type = c("n", "c"), return_order = TRUE)
  
  for(i in 1:ncol(sorted_table)){
    geno_for_pairscan <- geno_for_pairscan[,sorted_table[,i]]
  }
  
  
  if(verbose){cat("Checking for linear independence...\n")}
  data_obj$geno_for_pairscan <- geno_for_pairscan
  geno_ind <- get_linearly_independent(data_obj)
  
  rownames(geno_ind$independent_markers) <- rownames(data_obj$pheno)
  data_obj$geno_for_pairscan <- geno_ind$independent_markers
  data_obj$effect_size_cutoff <- min_effect_size
  data_obj$peak_density = peak_density
  data_obj$window_size = window_size
  data_obj$tolerance = tolerance
  
  if(verbose){
    cat(length(geno_ind[[2]]), "allele(s) rejected.\n")
    cat("Final alleles selected:", "\t", ncol(geno_ind$independent_markers), "\n")
  }

  if(plot_peaks){
    dev.off()
  }
  
  return(data_obj)
  
  
}