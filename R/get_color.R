#
#plotting. Given a color name it returns the hex 
#colors used to make color ramps

#' get a hex color string
#' 
#' This internal function stores colors for plotting.
#' Given a color name ("green", "purple", "red", "orange", 
#' "blue", "brown", "yellow", "gray")
#' and a darkness: "l" for light colors, "d" for dark colors 
#' and "f" for the full light/dark spectrum
#' this function returns the corresponding ramp of colors.
#'
#' @param col_name string color name. Must be one of "green", 
#' "purple", "red", "orange", "blue", "brown", "yellow", "gray"
#' @param light_dark character value. One of ("f", "l", "d").
#' "l" indicates light colors, "d" indicates dark colors, and "f"
#' indicates colors ranging from light to dark.
#'
#' @return a vector of length four containing the hex colors 
#' indicated by the parameters
#' 
#'
#' @keywords internal

get_color <- function(col_name = c("green", "purple", "red", "orange", "blue", "brown", "yellow", "gray"), 
light_dark = c("f", "l", "d")){
  
  
  col_name = col_name[1]
  
  light_dark_check <- grep("f", light_dark)
  if(length(light_dark_check) > 0){light_dark = "f"}
  
  possible_light_dark <- c("f", "l", "d")
  light_dark_check2 <- match(light_dark, possible_light_dark)
  if(is.na(light_dark_check2)){
    message("Possible specifications of light_dark are: ", paste(possible_light_dark, collapse = ","))
    stop()
  }
  
  possible_cols <- c("green", "purple", "red", "orange", "blue", "brown", "yellow", "gray")		
  col_check <- match(col_name, possible_cols)
  
  if(is.na(col_check)){
    message("Possible colors are: ", paste(possible_cols, collapse = ", "))
    stop()
  }
  
  
  light_mat <- matrix(
    c("#edf8fb", "#ccece6", "#99d8c9", "#66c2a4",
      "#f2f0f7", "#dadaeb", "#bcbddc", "#9e9ac8",
      "#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a",
      "#feedde", "#fdd0a2", "#fdae6b", "#fd8d3c",
      "#eff3ff", "#c6dbef", "#9ecae1", "#6baed6",
      "#f5f5f5", "#f6e8c3", "#dfc27d", "#bf812d",
      "#ffffe5", "#fff7bc", "#fee391", "#fec44f",
      "#ffffff", "#f0f0f0", "#d9d9d9", "#bdbdbd"), nrow = 4, byrow = FALSE)
  
  
  dark_mat <- matrix(
    c("#66c2a4", "#41ae76", "#238b45", "#005824", #green
      "#9e9ac8", "#807dba", "#6a51a3", "#4a1486", #purple
      "#fb6a4a", "#ef3b2c", "#cb181d", "#99000d", #red
      "#fd8d3c", "#f16913", "#d94801", "#8c2d04",
      "#6baed6", "#4292c6", "#2171b5", "#084594",
      "#dfc27d", "#bf812d", "#8c510a", "#543005",
      "#fec44f", "#fe9929", "#ec7014", "#cc4c02",
      "#bdbdbd", "#969696", "#737373", "#525252"), nrow = 4, byrow = FALSE)
  
  
  full_mat <- rbind(light_mat[c(1,3),], dark_mat[c(2,4),])
  
  if(light_dark == "l"){
    # color order: same as in arguments, light to dark 
    all_col_ref <- light_mat
  }
  if(light_dark == "d"){
    #color order: same as in arguments, light to dark 
    all_col_ref <- dark_mat
  }
  
  if(light_dark == "f"){
    all_col_ref <- 	full_mat
  }
  
  
  colnames(all_col_ref) <- possible_cols
  
  col_locale <- which(colnames(all_col_ref) == col_name)
  
  return(all_col_ref[,col_locale])
  
}