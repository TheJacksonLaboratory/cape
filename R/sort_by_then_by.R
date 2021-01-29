#' Sort a table by a list of columns
#' 
#' This internal function sorts a table first by one column,
#' and then by a second column.
#' Columns can contain numeric values or string values. These
#' are specified for each column by col_type.
#' The values of the second column are sorted by successive 
#' values in the first column.
#' if return order is TRUE, the function returns the order of
#' each column instead of the ordered table. To order the original
#' table, sort first by the order in the first column of the order
#' matrix and then by the offer in the second colum of the order
#' matrix.
#' 
#' @param tableX a matrix
#' @param sort_cols A vector of length two, indicating which
#' columns in the matrix to sort by, and in which order. To
#' sort first by column 1 and then by column 2, sort_cols should
#' be c(1,2).
#' @param col_type A vector of length two indicating whether each
#' column contains character values ("c") or numeric values ("n").
#' Specified in the same order as sort_cols.
#' @param decreasing Whether values should be sorted in decreasing order
#' @param return_order Whether to return the sorted table (FALSE) or to 
#' return the order used to sort the table (TRUE). Defaults to FALSE.
#' 
#' @return If return_order is FALSE, this function returns tableX sorted by sort_cols.
#' If return_order is TRUE, this function returns a two-column matrix indicating the order in which
#' to put tableX to sort it by sort_cols. To sort a table to match the order, order first by the
#' first column and then by the second column.
#' @keywords internal
#' 
sort_by_then_by <- function(tableX, sort_cols = c(1,2), col_type = c("c", "n"), decreasing = FALSE, return_order = FALSE){
  
  if(length(sort_cols) > 2){
    stop("This script can't handle more than 2 columns")
  }
  
  if(length(sort_cols) != length(col_type)){
    stop("The col_type vector must be the same length as sort_cols")
  }
  
  if(length(sort_cols) == 1){
    
    if(col_type == "n"){
      sorted_col <- sort.int(as.numeric(tableX[,sort_cols]), index.return = TRUE, decreasing = decreasing, na.last = TRUE, method = "radix")
    }else{
      sorted_col <- sort.int(tableX[,sort_cols], index.return = TRUE, decreasing = decreasing, na.last = TRUE, method = "radix")
    }
    
    new_table <- tableX[sorted_col$ix,]
    return(new_table)
    
  }else{
    
    #start with the table ordered by the first sort column. We will change
    #chunks of this as we go
    if(col_type[1] == "n"){
      new_order <- sort.int(as.numeric(tableX[,sort_cols[1]]), index.return = TRUE, na.last = TRUE, decreasing = decreasing, method = "radix")$ix	
    }else{
      new_order <- sort.int(tableX[,sort_cols[1]], index.return = TRUE, decreasing = decreasing, na.last = TRUE, method = "radix")$ix
    }
    sorted_table <- tableX[new_order, ]
    
    if(return_order){
      final_order <- new_order
    }
    
    if(col_type[1] == "n"){
      u_col_el <- sort(as.numeric(unique(sorted_table[,sort_cols[1]])), na.last = TRUE, method = "radix") #find the unique column elements for column 1
    }else{
      u_col_el <- sort(unique(sorted_table[,sort_cols[1]]), na.last = TRUE, method = "radix") #find the unique column elements for column 1
    }
    
    new_order <- NULL
    
    for(i in 1:length(u_col_el)){ #go through each of these elements, and sort the second column within the category of the first column
      if(is.na(u_col_el[i])){
        el_locale <- which(is.na(sorted_table[,sort_cols[1]])) #find the entries for element i in column 1
      }else{
        el_locale <- which(sorted_table[,sort_cols[1]] == u_col_el[i])
      }
      if(col_type[2] == "n"){
        subset_order <- sort.int(as.numeric(sorted_table[el_locale, sort_cols[2]]), index.return = TRUE, decreasing = decreasing, na.last = TRUE, method = "radix")$ix
      }else{
        subset_order <- sort.int(sorted_table[el_locale, sort_cols[2]], index.return = TRUE, decreasing = decreasing, na.last = TRUE, method = "radix")$ix
      }
      new_order <- c(new_order, el_locale[subset_order])
    }
    
    sorted_table <- sorted_table[new_order,]
    if(return_order){
      final_order <- cbind(final_order, new_order)
      return(final_order)
    }	
    
    return(sorted_table)
  } #end case for if sort_cols has more than one element

}
