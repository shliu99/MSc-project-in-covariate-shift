# function for supervised PCA

source("scripts/fun_kernel.R")


spca <- function(X,Y){
  X = as.matrix(X);
  X = t(X)
  n = ncol(X)
  L = linearkernel(Y)
  eig = eigen(L)
  delta = t(eig$vectors %*% diag(sqrt(ifelse(eig$values >= 0, eig$values, 0))))
  H = diag(n) - 1/n * outer(rep(1,n), rep(1,n))
  psi = X %*% H %*% t(delta)
  svd_res <- svd(psi)
  pc_n <- length(which(cumsum(svd_res$d)/sum(svd_res$d) <= 0.95)) + 1
  u <- psi %*% svd_res$v %*% diag(1 / svd_res$d)[,1:pc_n]
  return(u)
}
