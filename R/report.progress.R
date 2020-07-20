#' Report Progress of a Process
#' 
#' This function prints out the percent progress in a process
#' given the current iteration the total number of iterations
#' and the percentage at which a progress report is wanted.
#' It only works on loops in which the index is explicitly incremented
#' in the code.
#' 
#' @param current The current index being processed
#' @param total The total number of indices that will be processed
#' @param percent.text The percent progress at which the percent should
#' be printed to the screen.
#' @param percent.dot The percent progress at which a dot should be printed
#' to the screen
#' 
#' @return None. Prints output to the screen.


report.progress <- function(current, total, percent.text = 10, percent.dot = 2){
  
  all.iterations <- 1:total
  percent.prog <- round(all.iterations/total, 2)*100
  
  rounded.percent.write <- percent.prog%/%percent.text
  rounded.percent.dot <- percent.prog%/%percent.dot
  
  current.locale <- which(all.iterations == current)
  current.percent.write <- rounded.percent.write[current.locale]
  current.percent.dot <- rounded.percent.dot[current.locale]			
  
  curr.percent.write.locale <- which(rounded.percent.write == current.percent.write)
  curr.percent.dot.locale <- which(rounded.percent.dot == current.percent.dot)
  
  
  if(current.locale == min(curr.percent.write.locale)){
    cat(rounded.percent.write[current.locale]*percent.text, "%.", sep = "")
  }else{
    if(current.locale == min(curr.percent.dot.locale)){
      cat(".")
    }					
  }
  
}