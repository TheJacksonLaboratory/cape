#This is an internal function to find the
#column numbers of phenotypes when names
#are put in instead of column numbers
#if warn is set to TRUE, the default, the script
#prints a message if it can't find the column headers
#you are looking for. Setting warn to false is useful
#when you are using this script internally and just checking
#different matrices for column headers
get.col.num <- function (data.mat, col.which = NULL, dim.which = 2) 
{
  if (is.null(col.which)) {
    col.num <- 1:dim(data.mat)[dim.which]
  }
  else {
    if (is.numeric(col.which[1])) {
      return(col.which)
    }
    else {
      col.mat <- matrix(col.which, ncol = 1)
      find.col <- function(col.name, data.names) {
        col.num <- which(data.names == col.name)
        if(length(col.num) == 0){
          return(0)
        }
        return(col.num)
      }
      
      col.num <- unlist(apply(col.mat, 1, function(x) find.col(x, dimnames(data.mat)[[dim.which]])))
      didnt.find <- which(col.num == 0)
      if(length(didnt.find) > 0) {
        message("The following headers are not present:")
        cat(col.which[didnt.find], sep = "\n")
        col.num <- col.num[-didnt.find]
      }
    }
  }
  return(col.num)
}
