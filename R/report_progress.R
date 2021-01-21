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
#' @param percent_text The percent progress at which the percent should
#' be printed to the screen.
#' @param percent_dot The percent progress at which a dot should be printed
#' to the screen
#' @param verbose A logical value indicating whether to report progress to 
#' the screen. Defaults to TRUE.
#' 
#' @return None. Prints output to the screen.
#' @keywords internal


report_progress <- function(current, total, percent_text = 10, percent_dot = 2, verbose = TRUE){
  
  all_iterations <- 1:total
  percent_prog <- round(all_iterations/total, 2)*100
  
  rounded_percent_write <- percent_prog%/%percent_text
  rounded_percent_dot <- percent_prog%/%percent_dot
  
  current_locale <- which(all_iterations == current)
  current_percent_write <- rounded_percent_write[current_locale]
  current_percent_dot <- rounded_percent_dot[current_locale]			
  
  curr_percent_write_locale <- which(rounded_percent_write == current_percent_write)
  curr_percent_dot_locale <- which(rounded_percent_dot == current_percent_dot)
  
  
  if(current_locale == min(curr_percent_write_locale)){
    if(verbose){cat(rounded_percent_write[current_locale]*percent_text, "%.", sep = "")}
  }else{
    if(current_locale == min(curr_percent_dot_locale)){
      if(verbose){cat(".")}
    }
  }
  
}