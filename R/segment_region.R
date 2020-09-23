#' Divide a region into equal parts.
#' 
#' This is an internal function used to segment regions 
#' for plotting. It returns n evenly spaced points in a given
#' region. The points can be aligned to the ends of 
#' the interval, or centered in the interval.
#' 
#' @param region_min A numerical value indicating the minimum value of the region.
#' @param region_max A numerical value indicating the maximum value of the region.
#' @param num_points The number of points to place in the region.
#' @param alignment A character element indicating whether the points should be 
#' centered within the region or whether they should extend to the ends of the region.
#' 
#' @return Returns n points spaced evenly across the defined region
#'
segment_region <- function(region_min, region_max, num_points, alignment = c("center", "ends")){
  
  if(num_points < 2){
    return(mean(c(region_min, region_max)))
  }
  
  if(length(grep("c", alignment)) > 0){
    alignment <- "center"
  }
  
  
  total_region <- region_max - region_min
  
  if(alignment == "ends"){
    point_seq <- seq(region_min, region_max, total_region/(num_points-1))
    return(point_seq)
  }
  
  
  if(alignment == "center"){
    #first break the segment into n+1 regions
    point_seq <- seq(region_min, region_max, total_region/num_points)
    #find the center of each region
    cons_pairs <- consec_pairs(1:length(point_seq))
    center_points <- apply(cons_pairs, 1, function(x) mean(c(point_seq[x[1]], point_seq[x[2]])))
    return(center_points)
  }
  
  
  # plot(center.points, rep(1, length(center.points)), xlim = c(region_min, region_max), col = "red")
  # points(point_seq, rep(1.2, length(point_seq)), col = "blue")
  
  
}
