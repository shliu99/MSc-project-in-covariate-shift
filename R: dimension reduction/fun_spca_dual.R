# function for supervised PCA - improved efficiency on svd
library(rARPACK)
source("scripts/fun_kernel.R")


spca <- function(X,Y, prop = 0.95){
  if(!is.matrix(X)){
    X = as.matrix(X)
  }
  X = t(X)
  n = ncol(X)
  L = linearkernel(Y)
  eig = eigen(L)
  rm(L)
  delta = t(eig$vectors %*% diag(sqrt(ifelse(eig$values >= 0, eig$values, 0))))
  H = diag(n) - 1/n * outer(rep(1,n), rep(1,n))
  psi = X %*% H %*% t(delta)
  svd_sing <- svd(psi, nu = 0, nv = 0)$d
  pc_n <- length(which(cumsum(svd_sing^2)/sum(svd_sing^2) <= prop)) + 1
  v = svds(psi, k = pc_n, nu = 0, nv = pc_n)$v
  if(pc_n >= 2){
    u <- psi %*% v %*% diag(1 / svd_sing[1:pc_n])
  } else{
    u <- psi %*% v %*% (1 / svd_sing[1:pc_n])
  }
  return(u)
}

