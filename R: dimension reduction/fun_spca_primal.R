# function for supervised PCA - primal
library(rARPACK)
source("scripts/fun_kernel.R")


spca <- function(X,Y, prop = 0.95){
  if(!is.matrix(X)){
    X = as.matrix(X)
  }
  X = t(X)
  n = ncol(X)
  L = linearkernel(Y)
  H = diag(n) - 1/n * outer(rep(1,n), rep(1,n))
  Q = X %*% H %*% L %*% H %*% t(X)
  if(det(Q) == 0){
    Q = Q + diag(1e-6,nrow(X))
  }
  eig.vals <- eigen(Q, only.values = TRUE)$values
  pc_n <- length(which(cumsum(eig.vals^2)/sum(eig.vals^2) <= 0.95)) + 1
  beta <- eigs(Q, k = pc_n, which = "LM")$vectors
  return(beta)
}

