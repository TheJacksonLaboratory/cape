#' Get DO colors
#' 
#' This function returns the official DO allele colors.
#' as described here:
#' \url{http://www.csbio.unc.edu/CCstatus/index.py?run=AvailableLines.information}
#' 
#' @param color_scheme "DO/CC" or "other" The default is "DO/CC"
#' If set to "other", colors unrelated to the DO/CC are used. 
#' @param alleles an vector af alleles in the set "AJ", "B6", "129", 
#' "NOD", "NZO", "CAST", "PWK", "WSB". This argument can be used to 
#' retrieve colors only for a subset of alleles.
#' 
#' @return A table with three columns corresponding to strain names,
#' strain nicknames, and the assigned color values.
#' 
#' @keywords internal
#' 
get_allele_colors <- function(color_scheme = c("DO/CC", "other"), alleles = NULL){
    
    
     color_scheme = color_scheme[1]
    
    if(color_scheme == "DO/CC"){
      aj <- rgb(240/256 ,240/256 ,0/256)
      b6 <- rgb(128/256,128/256,128/256)
      c <- rgb(240/256, 128/256, 128/256)
      nod <- rgb(16/256, 16/256, 240/256)
      nzo <- rgb(0/256, 160/256, 240/256)
      cast <- rgb(0/256, 160/256, 0/256)
      pwk <- rgb(240/256, 0/256, 0/256)
      wsb <- rgb(144/256, 0/256, 224/256)
      
      strain_names <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")	
      strain_nicknames <- LETTERS[1:8]
      allele_colors <- c(aj, b6, c, nod, nzo, cast, pwk, wsb)	
    }else{
      allele_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00")
      strain_names <- LETTERS[1:8]
      strain_nicknames <- 	LETTERS[1:8]
    }
    
    allele_table <- cbind(strain_names, strain_nicknames, allele_colors)
    
    
    if(!is.null(alleles)){
      allele_idx <- match(alleles, strain_names)
      if(!is.na(allele_idx[1])){
        allele_table <- allele_table[allele_idx,,drop=FALSE]
      }else{
        allele_idx <- match(alleles, strain_nicknames)
        if(!is.na(allele_idx[1])){
          allele_table <- allele_table[allele_idx,,drop=FALSE]
        }
      }
    }
    
    return(allele_table)
    
  }
