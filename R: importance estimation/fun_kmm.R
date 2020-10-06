# function for KMM

library("quadprog")
source("scripts/fun_kernel.R"); source("scripts/fun_gaussian_cv.R")

kmm <- function(x_tr,x_te){
  x_tr = as.matrix(x_tr); x_te = as.matrix(x_te)
  n_tr = nrow(x_tr)
  n_te = nrow(x_te)
  sigma = kernelWidth_median(x_tr)
  # calculate kernel
  cat("computing kernel for training data... \n")
  Ktrtr = rbfkernel(x_tr, sigma = sigma)
  Ktrtr = 0.9*(Ktrtr + t(Ktrtr))
  Xtrtr = 
  if (det(Ktrtr) <= 0) {
    Ktrtr = Ktrtr + diag(1e-5,n_tr)
  }
  if (det(Ktrtr) <= 0) {
    Ktrtr = Ktrtr + diag(1e-4,n_tr)
  }
  if (det(Ktrtr) <= 0) {
    Ktrtr = Ktrtr + diag(1e-3,n_tr)
  }
  if (det(Ktrtr) <= 0) {
    Ktrtr = Ktrtr + diag(1e-2,n_tr)
  }
  while(det(Ktrtr) <= 0) {
    Ktrtr = Ktrtr + diag(1e-1,n_tr)
  }
  # compute kappa
  cat("computing kappa... \n")
  Ktrte = rbfkernel(X = x_tr, Y = x_te, sigma = sigma)
  kappa = rowSums(Ktrte) * n_tr/n_te
  # compute epsilon
  cat("computing epsilon... \n")
  eps = (sqrt(n_tr) - 1)/sqrt(n_tr)
  # constraints
  Amat <- cbind(rep(-1, n_tr), rep(1,n_tr), diag(1, n_tr), diag(-1, n_tr))
  bvec <- c(-n_tr*(eps + 1), max(-n_tr*(eps - 1), 0), rep(0, n_tr), 
            rep(-sqrt(n_tr) + 1, n_tr))
  # computing beta
  cat("computing beta... \n")
  beta = solve.QP(Ktrtr, kappa, Amat, bvec)$solution
  beta[beta < 0] = 0
  return(beta)
}

