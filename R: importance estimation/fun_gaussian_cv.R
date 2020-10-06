# This script contains two methods for choosing Gaussian kernel width

# by HSIC
kernelWidth_hsic <- function(X, Y){ 
  X = as.matrix(X);
  n = nrow(X)
  H = diag(n) - 1/n * outer(rep(1,n), rep(1,n))
  vals = seq(0.01,1,0.01)
  hsic = matrix(NA, nrow = length(vals), ncol = 2)
  hsic[,1] <- seq(0.01,1,0.01)
  L = linearkernel(Y)
  for(i in 1:length(vals)){
    sig = vals[i]
    K = rbfkernel(X, sigma = sig)
    hsic[i,2] = sum(diag(K %*% H %*% L %*% H))
  }
  out = hsic[which(hsic[,2] == max(hsic[,2])), 1][1]
  return(out)
}

# by median distance
kernelWidth_median <- function(X, Y = NULL){
  n <- nrow(X)
  if (is.null(Y)) {
    XtX <- tcrossprod(X)
    XX <- matrix(1, n) %*% diag(XtX)
    D <- XX - 2 * XtX + t(XX)
    sigma = median(D)
  }else{
    m <- nrow(Y)
    XX <- matrix(apply(X ^ 2, 1, sum), n, m)
    YY <- matrix(apply(Y ^ 2, 1, sum), n, m, byrow = TRUE)
    XY <- tcrossprod(X, Y)
    D <- XX - 2 * XY + YY
    sigma = median(D)
  }
  return(sigma)
}

