#'
#' @export
inverse <- function(a,b=diag(1,ncol(a))){
  a <- as.matrix(a)
  m <- cbind(a, b)
  if(nrow(a)!= ncol(a))
    stop("Cannot obtain Inverse. Matrix not square",call. = F)
  v <- seq_len(nrow(a))
  for(i in v){
    if(m[i,i] == 0)
      stop("rank less than ",ncol(a),call. = F)
    m[i,] <- m[i,]/m[i,i]
    for(j in v[-i]){
      m[j,] <-m[j,] -  m[i,]*m[j,i]
    }
  }
  m[,-v]
}

#' @export
mat_det <- function(m){
  if(length(m) == 1) return(m)
  n <- 0
  for(i in seq_along(m[,1])){
    n <- n + (-1)^(i-1) * m[1,i] * mat_det(m[-1,-i])
  }
  n
}

ei <- function(x, A){
  mat_det(A - diag(x,nrow(A)))
}











