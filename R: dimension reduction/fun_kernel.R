# This script contains two kernel functions - rbf and linear

# rbfkernel function - default sigma = median distance

rbfkernel <- function(X, sigma = NULL, Y = NULL) {
  if (!is.matrix(X))
  {print("X must be a matrix containing samples in its rows")
    return()
  }
  if (!is.null(Y)){
    if (!is.matrix(Y))
    {print("Y must be a matrix containing samples in its rows or NULL if it should not be used")
      return()
    }
    if (ncol(X) != ncol(Y)){
      print("The samples in the rows of X and Y must be of same dimension")
      return()
    }
  }
  n <- nrow(X)
  if (is.null(Y)) {
    XtX <- tcrossprod(X)
    XX <- matrix(1, n) %*% diag(XtX)
    D <- XX - 2 * XtX + t(XX)
    if(is.null(sigma)){
      sigma = var(D)
    }
  }
  else{
    m <- nrow(Y)
    XX <- matrix(apply(X ^ 2, 1, sum), n, m)
    YY <- matrix(apply(Y ^ 2, 1, sum), n, m, byrow = TRUE)
    XY <- tcrossprod(X, Y)
    D <- XX - 2 * XY + YY
    if(is.null(sigma)){
      sigma = var(D)
    }
  }
  # calculate rbf-kernel matrix
  K <- exp(-D / (2 * sigma))
  return(K)
}



#linear kernel
linearkernel <- function(X, Y = NULL){
  X = as.matrix(X)
  if(is.null(Y)){
    out <- tcrossprod(X,X)
  } else{
    out <- tcrossprod(X, Y)
  }
  return(out)
}
