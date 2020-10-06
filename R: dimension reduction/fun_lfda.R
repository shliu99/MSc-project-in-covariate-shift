# function for lfda supervised dimension reduction

library("lfda");

lfda <- function (x, y, r, metric = c("orthonormalized", "plain", 
                                      "weighted"), knn = 7){
  metric <- match.arg(metric)
  x <- t(as.matrix(x))
  y <- t(as.matrix(y))
  d <- nrow(x)
  n <- ncol(x)
  if (is.null(r)) 
    r <- d
  tSb <- mat.or.vec(d, d)
  tSw <- mat.or.vec(d, d)
  for (i in unique(as.vector(t(y)))) {
    Xc <- x[, y == i]
    nc <- ncol(Xc)
    Xc2 <- t(as.matrix(colSums(Xc^2)))
    distance2 <- repmat(Xc2, nc, 1) + repmat(t(Xc2), 1, nc) - 
      2 * t(Xc) %*% Xc
    A <- getAffinityMatrix(distance2, knn, nc)
    Xc1 <- as.matrix(rowSums(Xc))
    G <- Xc %*% (repmat(as.matrix(colSums(A)), 1, d) * t(Xc)) - 
      Xc %*% A %*% t(Xc)
    tSb <- tSb + (G/n) + Xc %*% t(Xc) * (1 - nc/n) + Xc1 %*% 
      (t(Xc1)/n)
    tSw <- tSw + G/nc
  }
  X1 <- as.matrix(rowSums(x))
  tSb <- tSb - X1 %*% t(X1)/n - tSw
  tSb <- (tSb + t(tSb))/2
  tSw <- (tSw + t(tSw))/2
  if(det(tSw) == 0){
    tSw = tSw + diag(1e-6,d)
  }
  if (r == d) {
    eigTmp <- eigen(solve(tSw) %*% tSb)
  }
  else {
    eigTmp <- suppressWarnings(rARPACK::eigs(A = solve(tSw) %*% 
                                               tSb, k = r, which = "LM"))
  }
  eigVec <- Re(eigTmp$vectors)
  eigVal <- as.matrix(Re(eigTmp$values))
  Tr <- getMetricOfType(metric, eigVec, eigVal, d)
  Z <- t(t(Tr) %*% x)
  out <- list(T = Tr, Z = Z)
  class(out) <- "lfda"
  return(out)
}
