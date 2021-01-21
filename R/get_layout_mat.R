#' Get the best layout matrix for a given number of panes per page.
#' 
#' This function is for automatically deciding how to lay out multiple
#' images per page. It takes as an argument the number of images the user
#' wants on a single page and returns the most square matrix possible.
#' 
#' @param num_panes The number of images to plot per page.
#' @param type A character vector specifying whether the layout matrix 
#'   should have more columns than rows ("landscape"), or more rows than 
#'   columns ("upright"), if the resulting matrix is not square.
#'
#' @return A layout matrix with positions for each element in num_panes
#' and zeros filling the rest of the rectangle where nothing will be plotted.
#' @keywords internal
#' 

get_layout_mat <- function(num_panes, type = c("landscape", "upright")){
    
    #first find the square root of the number of panes
    edge_dim <- sqrt(num_panes)
    round_dim <- round(edge_dim)
    
    #if the number is square, just return a square matrix
    if(round_dim == edge_dim){
      layout_mat <- matrix(1:num_panes, ncol = edge_dim, byrow = TRUE)
      return(layout_mat)
    }
    
    #if the the edge_dim is not a whole number
    # we need to do some more figuring
    #if the nearest whole number is less than
    #the square root of the number of panes
    #add one to get close
    if(round_dim < edge_dim){
      dims <- c(round_dim, (round_dim + 1))
    }else{
      dims <- c(round_dim, (round_dim - 1)) #otherwise, subtract one to get close
    }
    
    total_panes <- prod(dims)
    test_dims <- matrix(rep(dims, 2), nrow = 2, byrow = TRUE)
    
    #as long as we haven't found a spot for all plots
    while(all(total_panes < num_panes)){
      #add rows and columns until we get
      #equal to or greater than the total
      #number of panes.
      for(i in 1:length(test_dims[,1])){
        max_dim_locale <- which(test_dims[1,] == max(test_dims[1,]))
        test_dims[i,max_dim_locale[1]] <- test_dims[i,max_dim_locale[1]] + 1
      }		
      total_panes <- apply(test_dims, 1, prod)
    }
    
    
    #figure out if either of these two options equals the original num_panes
    found_it <- which(total_panes == num_panes)
    if(length(found_it) > 0){
      found_it_locale <- which(total_panes == num_panes)
      dims <- test_dims[found_it_locale[1],]
    }else{
      #otherwise, figure out which of these is the minimum of 
      #the values that goes over the total number of panes
      greater <- which(total_panes > num_panes)
      dim_locale <- which(total_panes == min(total_panes[greater]))
      dims <- test_dims[dim_locale[1],]
    }
    
    
    #Make the layout matrices
    #fill in 0's for the panes we don't
    #want to fill with figures
    
    #the default layout is landscape
    if(length(grep("l", type)) > 0){
      layout_mat <- matrix(c(1:num_panes, rep(0, (total_panes[1]-num_panes))), ncol = max(dims), byrow = TRUE)
    }else{
      layout_mat <- matrix(c(1:num_panes, rep(0, (total_panes[1]-num_panes))), ncol = min(dims), byrow = TRUE)
    }
    
    return(layout_mat)
  }

