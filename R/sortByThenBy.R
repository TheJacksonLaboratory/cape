#' Sort a table by a list of columns
#' 
#' This internal function sorts a table first by one column,
#' and then by a second column.
#' Columns can contain numeric values or string values. These
#' are specified for each column by col.type.
#' The values of the second column are sorted by successive 
#' values in the first column.
#' if return order is TRUE, the function returns the order of
#' each column instead of the ordered table. To order the original
#' table, sort first by the order in the first column of the order
#' matrix and then by the ofer in the second colum of the order
#' matrix.
#' 
#' @param tableX a matrix
#' @param sort.cols A vector of length two, indicating which
#' columns in the matrix to sort by, and in which order. To
#' sort first by column 1 and then by column 2, sort.cols should
#' be c(1,2).
#' @param col.type A vector of length two indicating whether each
#' column contains character values ("c") or numeric values ("n").
#' Specified in the same order as sort.cols.
#' @param decreasing Whether values should be sorted in decreasing order
#' @param return.order Whether to return the sorted table (FALSE) or to 
#' return the order used to sort the table (TRUE). Defaults to FALSE.
#' 
#' @return If return.order is FALSE, this function returns tableX sorted by sort.cols.
#' If return.order is TRUE, this function returns a two-column matrix indicating the order in which
#' to put tableX to sort it by sort.cols. To sort a table to match the order, order first by the
#' first column and then by the second column.
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
