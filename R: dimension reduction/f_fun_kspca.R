# function for kernel SPCA

source("scripts/fun_kernel.R")
library(rARPACK)

kSPCA <- function(X, Y, sigma){
  X = as.matrix(X); Y = as.matrix(Y)
  K = rbfkernel(X, sigma = sigma)
  L = linearkernel(Y)
  n = nrow(K)
  if(det(K) == 0){
    K = K + diag(1e-6, n)
  }
  H = diag(n) - 1/n * outer(rep(1,n), rep(1,n))
  Q = K %*% L %*% H %*% K
  
  # eig.res <- eigen(Q)
  eig.vals <- eig.res$values
  real.eigs.idx <- which(abs(Im(eig.vals)) < 1e-6)
  real.eigs <- abs(Re(eig.vals[real.eigs.idx]))
  ordered.real.eigs <- real.eigs[order(real.eigs, decreasing = T)]
  ord <- real.eigs.idx[order(real.eigs, decreasing = T)]
  pc_n <- length(which(cumsum(ordered.real.eigs)/sum(ordered.real.eigs) <= 0.95)) + 1
  selected.real.idx <- ord[1:pc_n]
  eig.vec <- Re(eigs(Q, max(selected.real.idx))[,selected.real.idx])
  # beta <- Re(eig.res$vectors[,selected.real.idx])
  return(beta)
}

