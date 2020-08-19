#' Bins a single scan curve into peaks automiatically
#' 
#' This is an internal function used to select markers
#' for the pair scan based on single scan results. The 
#' algorithm first finds the difference between all 
#  consecutive points in the vector. It then looks for 
#' runs of all positive and all negative values.
#' It smooths the curve and identifies peaks exceeding
#' the threshold defined by amp.min.
#' 
#' @param the.curve vector representing the curve to be binned into peaks
#' @param plot.peaks default = FALSE
#' @param window.size A numeric value setting how many markers 
#' should be included in each window. If NULL, the window size is
#' set to the maximum number of consecutive rises or falls in the 
#' curve
#' @param amp.min A numeric value indicating the minimum magnitude
#' of a peak. All peaks below this magnitude will be removed. If 
#' NULL amp.min is set to the sd of the curve/2.
#'
#' @return This function returns a list with the following elements:
#' bins: a vector the same length as the input curve identifying which
#' peak each position was assigned to. 
#' smoothed.curve: A vector defining the smoothed curve
#' window.size: The input window.size or the cacluated window.size if 
#' window.size was NULL
#' amp.min: the input amp.min or calculated amp.min if amp.min was NULL

bin.curve <- function(the.curve, plot.peaks = FALSE, window.size = NULL, amp.min = NULL){
  
    if(all(is.na(the.curve))){
      return(list("bins" = the.curve))
    }
  
  #====================================================================		
  # internal functions
  #====================================================================		
  #This function finds indices of two peaks between two troughs
  #or vice versa
  #The function returns indices in the index.by vector
  #that mark the beginning of a run of multiple elements
  #of the same state.
  remove.runs <- function(peak.locale, trough.locale, smoothed.curve){
    bigv <- c(peak.locale, trough.locale)
    bigc <- c(rep(1, length(peak.locale)), rep(2, length(trough.locale)))
    bigi <- c(1:length(peak.locale), 1:length(trough.locale))
    big.mat <- cbind(bigv, bigc, bigi)
    ind.order <- order(big.mat[,1])
    ordered.mat <- big.mat[ind.order,]
    
    consec.runs <- rle(ordered.mat[,2])
    run.locale <- which(consec.runs[[1]] > 1)
    
    remove.peaks <- NULL
    remove.troughs <- NULL
    
    #if there are runs, remove them and select the most
    #extreme point in each run to be the final point
    #high points for peaks, low points for troughs
    if(length(run.locale) > 0){
      run.idx <- vector(mode = "list", length = length(consec.runs[[1]]))
      start.num <- 1
      for(i in 1:length(run.idx)){
        run.idx[[i]] <- start.num:(start.num+consec.runs[[1]][i]-1)
        start.num <- max(unlist(run.idx))+1
      }
      
      for(i in 1:length(run.idx)){
        if(consec.runs[[1]][i] > 1){
          is.peak <- consec.runs[[2]][i] == 1
          curve.idx <- ordered.mat[run.idx[[i]],1]
          vals <- smoothed.curve[curve.idx]
          # plot(vals, type = "l")
          vector.idx <- ordered.mat[run.idx[[i]],3]
          
          if(is.peak){
            selected.idx <- which.max(vals)
            to.remove <- vector.idx[-selected.idx]
            remove.peaks <- c(remove.peaks, to.remove)
          }else{
            selected.idx <- which.min(vals)	
            to.remove <- vector.idx[-selected.idx]
            remove.troughs <- c(remove.troughs, to.remove)
          }
        }
      }
    }
    if(length(remove.peaks) > 0){
      peak.locale <- peak.locale[-remove.peaks]
    }
    if(length(remove.troughs) > 0){
      trough.locale <- trough.locale[-remove.troughs]
    }
    
    trimmed.vectors <- list("peak.locale" = peak.locale, "trough.locale" = trough.locale)
    return(trimmed.vectors)
  }
  
  #====================================================================
  
  
  ymax = max(abs(the.curve), na.rm = TRUE)
  
  if(is.null(window.size)){
    all.diff <- diff(the.curve)
    pos.slope <- which(all.diff > 0)
    neg.slope <- which(all.diff < 0)
    pos.slope.runs <- diff(pos.slope)
    neg.slope.runs <- diff(neg.slope)
    max.pos <- max(pos.slope.runs)
    max.neg <- max(neg.slope.runs)
    window.size = min(c(max.pos, max.neg))
  }
  
  
  cols <- c("grey", "white")
  the.curve <- abs(the.curve)
  smoothed.curve <- caTools::runmean(the.curve, window.size)
  curve.bins <- rep(NA, length(the.curve))
  
  # hist(smoothed.curve)	
  if(is.null(amp.min)){
    amp.min = sd(smoothed.curve)/2
  }
  
  if(is.null(ymax)){ymax = max(abs(the.curve), na.rm = TRUE)}
  smoothed.x <- 1:length(smoothed.curve)
  
  smoothed.y <- c(0, smoothed.curve, 0)
  
  smoothed.left <- smoothed.curve - head(smoothed.y, length(smoothed.curve))
  smoothed.right <- smoothed.curve - tail(smoothed.y, length(smoothed.curve))
  
  peak.locale <- intersect(which(smoothed.left > 0), which(smoothed.right > 0))
  trough.locale <- intersect(which(smoothed.left < 0), which(smoothed.right < 0))
  
  #remove any runs of peaks and troughs, which happens when there are flat spots
  #in the smoothed curve
  trimmed.locale <- remove.runs(peak.locale, trough.locale, smoothed.curve)
  peak.locale <- trimmed.locale$peak.locale
  trough.locale <- trimmed.locale$trough.locale
  
  peak.y <- smoothed.curve[peak.locale]	
  
  trough.y <- smoothed.curve[trough.locale]	
  
  if(plot.peaks){
    par(mfrow = c(3,1))		
    # start.pt <- 400;end.pt = 500
    start.pt <- 1; end.pt <- length(smoothed.x)
    plot(x = start.pt:end.pt, the.curve[start.pt:end.pt], type = "l")
    points(smoothed.x[start.pt:end.pt], smoothed.curve[start.pt:end.pt], type = "l", col = "purple")
    points(smoothed.x[peak.locale], smoothed.curve[peak.locale], pch = 16, col = "red")
    points(smoothed.x[trough.locale], smoothed.curve[trough.locale], pch = 16, col = "blue")
  }
  
  
  #delete peaks and troughs that do not have high
  #enough amplitude
  #first look at peak to trough distances
  if(length(trough.y) < length(peak.y)){
    peak.trough.dist <- peak.y - c(trough.y, 0)
  }else{
    peak.trough.dist <- peak.y - trough.y	
  }
  
  small.drop <- which(peak.trough.dist < amp.min)
  #delete the troughs for each small drop and
  #delete any resulting runs in peaks
  if(length(small.drop) > 0){
    trough.locale <- trough.locale[-small.drop]
    if(length(trough.locale) > 0){
      trimmed.locale <- remove.runs(peak.locale, trough.locale, smoothed.curve)
      peak.locale <- trimmed.locale$peak.locale
      trough.locale <- trimmed.locale$trough.locale
    }
  }
  
  if(plot.peaks){		
    plot(smoothed.x[start.pt:end.pt], smoothed.curve[start.pt:end.pt], type = "l")
    points(smoothed.x[peak.locale], smoothed.curve[peak.locale], pch = 16, col = "red")
    points(smoothed.x[trough.locale], smoothed.curve[trough.locale], pch = 16, col = "blue")
  }
  
  padded.trough <- unique(c(1, trough.locale, length(smoothed.curve)))
  
  if(plot.peaks){
    plot.new()
    plot.window(xlim = c(1,length(smoothed.x)), ylim = c(min(the.curve, na.rm = TRUE), max(the.curve, na.rm = TRUE)))
  }
  
  for(j in 1:(length(padded.trough)-1)){
    poly.ind <- padded.trough[j]:padded.trough[j+1]
    poly.x <- floor(smoothed.x[poly.ind[1]]):ceiling(smoothed.x[tail(poly.ind, 1)])
    poly.y <- the.curve[poly.x]
    poly.y <- c(0,poly.y, 0)
    poly.x <- c(poly.x[1], poly.x, tail(poly.x, 1))
    curve.bins[poly.x] <- rep(j, length(poly.x))
    if(plot.peaks){
      polygon(x = poly.x, y = poly.y, col = cols[j%%length(cols)])
    }
  }#end looping through bins
  if(plot.peaks){
    points(smoothed.x, smoothed.curve, type = "l", ylim = c(min(the.curve), ymax), col = "purple")
    points(smoothed.x[peak.locale], smoothed.curve[peak.locale], pch = 16, col = "red")
    points(smoothed.x[trough.locale], smoothed.curve[trough.locale], pch = 16, col = "blue")
    axis(2)
  }
  
  return(list("bins" = curve.bins, "smoothed.curve" = smoothed.y, "window.size" = window.size, "amp.min" = amp.min))
}
