# function for kernel SPCA

source("scripts/fun_kernel.R"); library("Matrix")

kSPCA <- function(X, Y, sigma){
  X = as.matrix(X); Y = as.matrix(Y)
  K = rbfkernel(X, sigma = sigma)
  L = linearkernel(Y)
  n = nrow(K)
  H = diag(n) - 1/n * outer(rep(1,n), rep(1,n))
  Q = K %*% L %*% H %*% K
  if(det(Q) == 0){
    Q = nearPD(Q)$mat
  }
  decom <- eigen(Q)
  pc_n <- length(which(cumsum(decom$d)/sum(decom$d) <= 0.95)) + 1
  beta <- as.matrix(decom$vectors[, 1:pc_n])
  return(beta)
}
