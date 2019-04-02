# rank Z noramilze
rz.transform <- function (x, jitter = FALSE){
  x = rank(x, na.last = "keep", ties.method = "average") / (length(x) + 1)
  return(qnorm(x))
}
