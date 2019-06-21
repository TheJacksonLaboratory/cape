#' Sort a table by a list of columns
#' 
#' This script takes in a table and sorts by columns with each
#' successive column in the sort.cols vector sorted within the
#' factors of the previous column. Only the final column can 
#' contain continuous variables. All others should be discrete.
#' col.type defines whether each column contains a number ("n")
#' or a character ("c")
#' if return order is TRUE, the function returns the order of
#' each column instead of the ordered table
#' 
#' @param tableX a table object
#' @param sort.cols integer, number of columns to sort on (can be either 1 or 2)
#' @param col.type options are "c" (character) or "n" (number)
#' @param decreasing default = FALSE
#' @param retrun.order defualt = FALSE
#' 
sortByThenBy <- function(tableX, sort.cols = c(1,2), col.type = c("c", "n"), decreasing = FALSE, return.order = FALSE){
  
  if(length(sort.cols) > 2){
    stop("This script can't handle more than 2 columns")
  }
  
  if(length(sort.cols) != length(col.type)){
    stop("The col.type vector must be the same length as sort.cols")
  }
  
  if(length(sort.cols) == 1){
    
    if(col.type == "n"){
      sorted.col <- sort.int(as.numeric(tableX[,sort.cols]), index.return = TRUE, decreasing = decreasing, na.last = TRUE, method = "radix")
    }else{
      sorted.col <- sort.int(tableX[,sort.cols], index.return = TRUE, decreasing = decreasing, na.last = TRUE, method = "radix")
    }
    
    new.table <- tableX[sorted.col$ix,]
    return(new.table)
    
  }else{
    
    #start with the table ordered by the first sort column. We will change
    #chunks of this as we go
    if(col.type[1] == "n"){
      new.order <- sort.int(as.numeric(tableX[,sort.cols[1]]), index.return = TRUE, na.last = TRUE, decreasing = decreasing, method = "radix")$ix	
    }else{
      new.order <- sort.int(tableX[,sort.cols[1]], index.return = TRUE, decreasing = decreasing, na.last = TRUE, method = "radix")$ix
    }
    sorted.table <- tableX[new.order, ]
    
    if(return.order){
      final.order <- new.order
    }
    
    if(col.type[1] == "n"){
      u_col_el <- sort(as.numeric(unique(sorted.table[,sort.cols[1]])), na.last = TRUE, method = "radix") #find the unique column elements for column 1
    }else{
      u_col_el <- sort(unique(sorted.table[,sort.cols[1]]), na.last = TRUE, method = "radix") #find the unique column elements for column 1
    }
    
    
    new.order <- NULL
    
    for(i in 1:length(u_col_el)){ #go through each of these elements, and sort the second column within the category of the first column
      if(is.na(u_col_el[i])){
        el.locale <- which(is.na(sorted.table[,sort.cols[1]])) #find the entries for element i in column 1
      }else{
        el.locale <- which(sorted.table[,sort.cols[1]] == u_col_el[i])
      }
      if(col.type[2] == "n"){
        subset.order <- sort.int(as.numeric(sorted.table[el.locale, sort.cols[2]]), index.return = TRUE, decreasing = decreasing, na.last = TRUE, method = "radix")$ix
      }else{
        subset.order <- sort.int(sorted.table[el.locale, sort.cols[2]], index.return = TRUE, decreasing = decreasing, na.last = TRUE, method = "radix")$ix
      }
      new.order <- c(new.order, el.locale[subset.order])
    }
    
    sorted.table <- sorted.table[new.order,]
    if(return.order){
      final.order <- cbind(final.order, new.order)
      return(final.order)
    }	
    
    return(sorted.table)
  } #end case for if sort.cols has more than one element
  
  
  
  
  
  
}
