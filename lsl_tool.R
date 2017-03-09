.vec <- function(x) {
  x <- as.matrix(as.vector(x))
  return(x)
}

.vech <- function(x) {
  x <- as.matrix(x[lower.tri(x, diag = T)])
  return(x)
}

.ltri <- function(x) {
  x[upper.tri(x == 1)] = 0
  return(x)
}

.is_est <- function(x){
  x <- is.na(x) | x == 1
  return(x)
} # 判定是否需要估計，如果X為1為需要估計，X為NA為懲罰項，均返回"TRUE"。X為0為不需要估計，返回"FALSE"。

.is_one <- function(x){
  x <- !is.na(x) & x == 1
  return(x)
}

