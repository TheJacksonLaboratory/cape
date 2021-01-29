#' Bins a single scan curve into peaks automatically
#' 
#' This is an internal function used to select markers
#' for the pair scan based on single scan results. The 
#' algorithm first finds the difference between all 
#  consecutive points in the vector. It then looks for 
#' runs of all positive and all negative values.
#' It smooths the curve and identifies peaks exceeding
#' the threshold defined by amp_min.
#' 
#' @param the_curve vector representing the curve to be binned into peaks
#' @param plot_peaks default = FALSE
#' @param window_size A numeric value setting how many markers 
#' should be included in each window. If NULL, the window size is
#' set to the maximum number of consecutive rises or falls in the 
#' curve
#' @param amp_min A numeric value indicating the minimum magnitude
#' of a peak. All peaks below this magnitude will be removed. If 
#' NULL amp_min is set to the sd of the curve/2.
#'
#' @return This function returns a list with the following elements:
#' bins: a vector the same length as the input curve identifying which
#' peak each position was assigned to. 
#' smoothed_curve: A vector defining the smoothed curve
#' window_size: The input window_size or the calculated window_size if 
#' window_size was NULL
#' amp_min: the input amp_min or calculated amp_min if amp_min was NULL
#' 
#' @import caTools
#' @importFrom stats sd
#' @importFrom utils head tail
#' @importFrom graphics axis par plot plot.new plot.window points polygon
#' @keywords internal

bin_curve <- function(the_curve, plot_peaks = FALSE, window_size = NULL, amp_min = NULL){
  
    if(all(is.na(the_curve))){
      return(list("bins" = the_curve))
    }
  
  #====================================================================		
  # internal functions
  #====================================================================		
  #This function finds indices of two peaks between two troughs
  #or vice versa
  #The function returns indices in the index.by vector
  #that mark the beginning of a run of multiple elements
  #of the same state.
  remove_runs <- function(peak_locale, trough_locale, smoothed_curve){
    if(length(peak_locale) == 0 || length(trough_locale) == 0){
      trimmed_vectors <- list("peak_locale" = peak_locale, "trough_locale" = trough_locale)
      return(trimmed_vectors)
    }
    bigv <- c(peak_locale, trough_locale)
    bigc <- c(rep(1, length(peak_locale)), rep(2, length(trough_locale)))
    bigi <- c(1:length(peak_locale), 1:length(trough_locale))
    big_mat <- cbind(bigv, bigc, bigi)
    ind_order <- order(big_mat[,1])
    ordered_mat <- big_mat[ind_order,]
    
    consec_runs <- rle(ordered_mat[,2])
    run_locale <- which(consec_runs[[1]] > 1)
    
    remove_peaks <- NULL
    remove_troughs <- NULL
    
    #if there are runs, remove them and select the most
    #extreme point in each run to be the final point
    #high points for peaks, low points for troughs
    if(length(run_locale) > 0){
      run_idx <- vector(mode = "list", length = length(consec_runs[[1]]))
      start_num <- 1
      for(i in 1:length(run_idx)){
        run_idx[[i]] <- start_num:(start_num+consec_runs[[1]][i]-1)
        start_num <- max(unlist(run_idx))+1
      }
      
      for(i in 1:length(run_idx)){
        if(consec_runs[[1]][i] > 1){
          is_peak <- consec_runs[[2]][i] == 1
          curve_idx <- ordered_mat[run_idx[[i]],1]
          vals <- smoothed_curve[curve_idx]
          # plot(vals, type = "l")
          vector.idx <- ordered_mat[run_idx[[i]],3]
          
          if(is_peak){
            selected_idx <- which.max(vals)
            to_remove <- vector.idx[-selected_idx]
            remove_peaks <- c(remove_peaks, to_remove)
          }else{
            selected_idx <- which.min(vals)	
            to_remove <- vector.idx[-selected_idx]
            remove_troughs <- c(remove_troughs, to_remove)
          }
        }
      }
    }
    if(length(remove_peaks) > 0){
      peak_locale <- peak_locale[-remove_peaks]
    }
    if(length(remove_troughs) > 0){
      trough_locale <- trough_locale[-remove_troughs]
    }
    
    trimmed_vectors <- list("peak_locale" = peak_locale, "trough_locale" = trough_locale)
    return(trimmed_vectors)
  }
  
  #====================================================================
  
  
  ymax = max(abs(the_curve), na.rm = TRUE)
  
  if(is.null(window_size)){
    all_diff <- diff(the_curve)
    pos_slope <- which(all_diff > 0)
    neg_slope <- which(all_diff < 0)
    pos_slope_runs <- diff(pos_slope)
    neg_slope_runs <- diff(neg_slope)
    if(length(pos_slope_runs) == 0 || length(neg_slope_runs) == 0){
      window_size = 1
    }else{
      max_pos <- max(pos_slope_runs)
      max_neg <- max(neg_slope_runs)
      window_size = min(c(max_pos, max_neg))
    }
  }
  
  
  cols <- c("grey", "white")
  the_curve <- abs(the_curve)
  smoothed_curve <- runmean(the_curve, window_size)
  curve_bins <- rep(NA, length(the_curve))
  
  # hist(smoothed_curve)	
  if(is.null(amp_min)){
    amp_min = sd(smoothed_curve)/2
  }
  
  if(is.null(ymax)){ymax = max(abs(the_curve), na.rm = TRUE)}
  smoothed_x <- 1:length(smoothed_curve)
  
  smoothed_y <- c(0, smoothed_curve, 0)
  
  smoothed_left <- smoothed_curve - head(smoothed_y, length(smoothed_curve))
  smoothed_right <- smoothed_curve - tail(smoothed_y, length(smoothed_curve))
  
  peak_locale <- intersect(which(smoothed_left > 0), which(smoothed_right > 0))
  trough_locale <- intersect(which(smoothed_left < 0), which(smoothed_right < 0))
  
  #remove any runs of peaks and troughs, which happens when there are flat spots
  #in the smoothed curve
  trimmed_locale <- remove_runs(peak_locale, trough_locale, smoothed_curve)
  peak_locale <- trimmed_locale$peak_locale
  trough_locale <- trimmed_locale$trough_locale
  
  peak_y <- smoothed_curve[peak_locale]	
  
  trough_y <- smoothed_curve[trough_locale]	
  
  if(plot_peaks){
    par(mfrow = c(3,1))		
    # start_pt <- 400;end_pt = 500
    start_pt <- 1; end_pt <- length(smoothed_x)
    plot(x = start_pt:end_pt, the_curve[start_pt:end_pt], type = "l")
    points(smoothed_x[start_pt:end_pt], smoothed_curve[start_pt:end_pt], type = "l", col = "purple")
    points(smoothed_x[peak_locale], smoothed_curve[peak_locale], pch = 16, col = "red")
    points(smoothed_x[trough_locale], smoothed_curve[trough_locale], pch = 16, col = "blue")
  }
  
  
  #delete peaks and troughs that do not have high
  #enough amplitude
  #first look at peak to trough distances
  if(length(trough_y) < length(peak_y)){
    peak_trough_dist <- peak_y - c(trough_y, 0)
  }else{
    peak_trough_dist <- peak_y - trough_y	
  }
  
  small_drop <- which(peak_trough_dist < amp_min)
  #delete the troughs for each small drop and
  #delete any resulting runs in peaks
  if(length(small_drop) > 0){
    trough_locale <- trough_locale[-small_drop]
    if(length(trough_locale) > 0){
      trimmed_locale <- remove_runs(peak_locale, trough_locale, smoothed_curve)
      peak_locale <- trimmed_locale$peak_locale
      trough_locale <- trimmed_locale$trough_locale
    }
  }
  
  if(plot_peaks){		
    plot(smoothed_x[start_pt:end_pt], smoothed_curve[start_pt:end_pt], type = "l")
    points(smoothed_x[peak_locale], smoothed_curve[peak_locale], pch = 16, col = "red")
    points(smoothed_x[trough_locale], smoothed_curve[trough_locale], pch = 16, col = "blue")
  }
  
  padded_trough <- unique(c(1, trough_locale, length(smoothed_curve)))
  
  if(plot_peaks){
    plot.new()
    plot.window(xlim = c(1,length(smoothed_x)), ylim = c(min(the_curve, na.rm = TRUE), max(the_curve, na.rm = TRUE)))
  }
  
  for(j in 1:(length(padded_trough)-1)){
    poly_ind <- padded_trough[j]:padded_trough[j+1]
    poly_x <- floor(smoothed_x[poly_ind[1]]):ceiling(smoothed_x[tail(poly_ind, 1)])
    poly_y <- the_curve[poly_x]
    poly_y <- c(0,poly_y, 0)
    poly_x <- c(poly_x[1], poly_x, tail(poly_x, 1))
    curve_bins[poly_x] <- rep(j, length(poly_x))
    if(plot_peaks){
      polygon(x = poly_x, y = poly_y, col = cols[j%%length(cols)])
    }
  }#end looping through bins
  if(plot_peaks){
    points(smoothed_x, smoothed_curve, type = "l", ylim = c(min(the_curve), ymax), col = "purple")
    points(smoothed_x[peak_locale], smoothed_curve[peak_locale], pch = 16, col = "red")
    points(smoothed_x[trough_locale], smoothed_curve[trough_locale], pch = 16, col = "blue")
    axis(2)
  }
  
  return(list("bins" = curve_bins, "smoothed_curve" = smoothed_y, "window_size" = window_size, "amp_min" = amp_min))
}
