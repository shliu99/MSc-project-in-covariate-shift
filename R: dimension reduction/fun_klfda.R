# function for kernel lfda supervised dimension reduction

library("lfda")

klfda <- function (k, y, r, metric = c("weighted", "orthonormalized", 
                                       "plain"), knn = 7, reg = 0.001){
  metric <- match.arg(metric)
  y <- t(as.matrix(y))
  n <- nrow(k)
  if (is.null(r)) 
    r <- n
  tSb <- mat.or.vec(n, n)
  tSw <- mat.or.vec(n, n)
  for (i in unique(as.vector(t(y)))) {
    Kcc <- k[y == i, y == i]
    Kc <- k[, y == i]
    nc <- nrow(Kcc)
    Kccdiag <- diag(Kcc)
    distance2 <- repmat(Kccdiag, 1, nc) + repmat(t(Kccdiag), 
                                                 nc, 1) - 2 * Kcc
    A <- getAffinityMatrix(distance2, knn, nc)
    Kc1 <- as.matrix(rowSums(Kc))
    Z <- Kc %*% (repmat(as.matrix(colSums(A)), 1, n) * t(Kc)) - 
      Kc %*% A %*% t(Kc)
    tSb <- tSb + (Z/n) + Kc %*% t(Kc) * (1 - nc/n) + Kc1 %*% 
      (t(Kc1)/n)
    tSw <- tSw + Z/nc
  }
  K1 <- as.matrix(rowSums(k))
  tSb <- tSb - K1 %*% t(K1)/n - tSw
  tSb <- (tSb + t(tSb))/2
  tSw <- (tSw + t(tSw))/2
  if(det(tSw) == 0){
    tSw = tSw + diag(1e-6, n)
  }
  if(det(tSb) == 0){
    tSb = tSb + diag(1e-6, n)
  }
  eigTmp <- suppressWarnings(rARPACK::eigs(A = solve(tSw + 
                                                       reg * diag(1, nrow(tSw), ncol(tSw))) %*% tSb, k = r, 
                                           which = "LM"))
  eigVec <- Re(eigTmp$vectors)
  eigVal <- as.matrix(Re(eigTmp$values))
  Tr <- getMetricOfType(metric, eigVec, eigVal, n)
  Z <- t(t(Tr) %*% k)
  out <- list(T = Tr, Z = Z)
  class(out) <- "lfda"
  return(out)
}
