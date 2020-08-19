#' Get DO colors
#' 
#' This function returns the official DO allele colors.
#' as described here \url{https://compgen.unc.edu/wp/?page_id=577}
#' 
#' @param color.scheme "DO/CC" or "other" The default is "DO/CC"
#' If set to "other", colors unrelated to the DO/CC are used. 
#' @param alleles an vector af alleles in the set "AJ", "B6", "129", 
#' "NOD", "NZO", "CAST", "PWK", "WSB". This argument can be used to 
#' retrieve colors only for a subset of alleles.
#' 
#' @return A table with three columns corresponding to strain names,
#' strain nicknames, and the assigned color values.
#' @export
get.allele.colors <-
  function(color.scheme = c("DO/CC", "other"), alleles = NULL){
    
    
     color.scheme = color.scheme[1]
    
    if(color.scheme == "DO/CC"){
      aj <- rgb(240/256 ,240/256 ,0/256)
      b6 <- rgb(128/256,128/256,128/256)
      c <- rgb(240/256, 128/256, 128/256)
      nod <- rgb(16/256, 16/256, 240/256)
      nzo <- rgb(0/256, 160/256, 240/256)
      cast <- rgb(0/256, 160/256, 0/256)
      pwk <- rgb(240/256, 0/256, 0/256)
      wsb <- rgb(144/256, 0/256, 224/256)
      
      strain.names <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")	
      strain.nicknames <- LETTERS[1:8]
      allele.colors <- c(aj, b6, c, nod, nzo, cast, pwk, wsb)	
    }else{
      allele.colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00")
      strain.names <- LETTERS[1:8]
      strain.nicknames <- 	LETTERS[1:8]
    }
    
    allele.table <- cbind(strain.names, strain.nicknames, allele.colors)
    
    
    if(!is.null(alleles)){
      allele.idx <- match(alleles, strain.names)
      if(!is.na(allele.idx[1])){
        allele.table <- allele.table[allele.idx,,drop=FALSE]
      }else{
        allele.idx <- match(alleles, strain.nicknames)
        if(!is.na(allele.idx[1])){
          allele.table <- allele.table[allele.idx,,drop=FALSE]
        }
      }
    }
    
    return(allele.table)
    
  }
