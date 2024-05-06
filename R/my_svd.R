

my_SVD <- function(A){
  n <- dim(A)
  transposed = FALSE
  if (n[1] < n[2]) {
    A <- t(A)
    transposed = TRUE
  }
  n <- dim(A)
  V <- matrix(0, n[2], n[2])
  U <- matrix(0, n[1], n[2])
  S <- numeric(n[2])
  for (i in seq.int(n[2])) {
    D <- crossprod(A)
    v <- D[1,]
    for(j in 1:20)  v <- D %*% v/v[1]
    V[,i] <- v / norm(v, "2")
    U[,i] <- A %*% V[,i]
    S[i] <- norm(U[,i], "2")
    U[,i] <- U[,i] / S[i]
    A <- A - S[i] * tcrossprod(U[,i], V[,i])
  }
  u <- U
  v <- V
  if(transposed){
    u <- V
    v <- U
  }
  list(d=S,u=u,v=v)
}
