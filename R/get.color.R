#This internal function stores colors for image 
#plotting. Given a color name it returns the hex 
#colors used to make color ramps


#' get a hex color string
#' 
#' given a color name ("green", "purple", "red", "orange", "blue", "brown", "yellow", "gray")
#' and a darkness ("f", "l", "d"), return the coressponding hex color string.
#'
#' @param col.name string color name
#' @param light.dar character in ("f", "l", "d")
#'
#' @return hex color string
#'
#' @export
get.color <- function(col.name = c("green", "purple", "red", "orange", "blue", "brown", "yellow", "gray"), light.dark = c("f", "l", "d")){
  
  
  col.name = col.name[1]
  
  light.dark.check <- grep("f", light.dark)
  if(length(light.dark.check) > 0){light.dark = "f"}
  
  possible.light.dark <- c("f", "l", "d")
  light.dark.check2 <- match(light.dark, possible.light.dark)
  if(is.na(light.dark.check2)){
    cat("Possible specifications of light.dark are:", possible.light.dark, sep = "\n")
    stop()
  }
  
  possible.cols <- c("green", "purple", "red", "orange", "blue", "brown", "yellow", "gray")		
  col.check <- match(col.name, possible.cols)
  
  if(is.na(col.check)){
    cat("Possible colors are:", possible.cols, sep = "\n")
    stop()
  }
  
  
  light.mat <- matrix(
    c("#edf8fb", "#ccece6", "#99d8c9", "#66c2a4",
      "#f2f0f7", "#dadaeb", "#bcbddc", "#9e9ac8",
      "#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a",
      "#feedde", "#fdd0a2", "#fdae6b", "#fd8d3c",
      "#eff3ff", "#c6dbef", "#9ecae1", "#6baed6",
      "#f5f5f5", "#f6e8c3", "#dfc27d", "#bf812d",
      "#ffffe5", "#fff7bc", "#fee391", "#fec44f",
      "#ffffff", "#f0f0f0", "#d9d9d9", "#bdbdbd"), nrow = 4, byrow = FALSE)
  
  
  dark.mat <- matrix(
    c("#66c2a4", "#41ae76", "#238b45", "#005824", #green
      "#9e9ac8", "#807dba", "#6a51a3", "#4a1486", #purple
      "#fb6a4a", "#ef3b2c", "#cb181d", "#99000d", #red
      "#fd8d3c", "#f16913", "#d94801", "#8c2d04",
      "#6baed6", "#4292c6", "#2171b5", "#084594",
      "#dfc27d", "#bf812d", "#8c510a", "#543005",
      "#fec44f", "#fe9929", "#ec7014", "#cc4c02",
      "#bdbdbd", "#969696", "#737373", "#525252"), nrow = 4, byrow = FALSE)
  
  
  full.mat <- rbind(light.mat[c(1,3),], dark.mat[c(2,4),])
  
  if(light.dark == "l"){
    # color order: same as in arguments, light to dark 
    all.col.ref <- light.mat
  }
  if(light.dark == "d"){
    #color order: same as in arguments, light to dark 
    all.col.ref <- dark.mat
  }
  
  if(light.dark == "f"){
    all.col.ref <- 	full.mat
  }
  
  
  colnames(all.col.ref) <- possible.cols
  
  col.locale <- which(colnames(all.col.ref) == col.name)
  
  return(all.col.ref[,col.locale])
  
}