# functions for implementing kernel mean match
# method described in paper: 
# http://is.tuebingen.mpg.de/fileadmin/user_upload/files/publications/shift-book-for-LeEtAl-webversion_5376[0].pdf
# section 1.3.3

library("quadprog")

matrix_dist <- function(mat1, mat2){
  k21 <- apply(mat2,1,function(x){(x-t(mat1))^2})
  ary <- array(k21, c(ncol(mat1), nrow(mat1), nrow(mat2)))
  out <- apply(ary, 2:3, sum)
  return(sqrt(out))
}

rbf_kernel <- function(mat1, mat2 = NULL,gamma = NULL){
  if(is.null(gamma)) {gamma = 1/ncol(mat1)}
  if(is.null(mat2)){
    K = dist(x, method = "euclidean", diag = TRUE, upper = TRUE)
  } else {
    K = matrix_dist(mat1, mat2)
  }
  out = exp(-gamma * K)
  return(out)
}


kmm <- function(x_tr,x_te,gamma = NULL){
  n_tr = nrow(x_tr)
  n_te = nrow(x_te)
  # calculate kernel
  cat("computing kernel for training data... \n")
  Ktrtr = rbf_kernel(x_tr, x_tr, gamma)
  # compute kappa
  cat("computing kappa... \n")
  Ktrte = rbf_kernel(x_tr, x_te, gamma)
  kappa = rowSums(Ktrte) * n_tr/n_te
  # compute epsilon
  cat("computing epsilon... \n")
  eps = (sqrt(n_tr) - 1)/sqrt(n_tr)
  # constraints
  Amat <- t(rbind(rep(-1, n_tr), rep(1,n_tr), diag(1, n_tr), diag(-1, n_tr)))
  bvec <- c(-n_tr*(eps + 1), -n_tr*(eps - 1), rep(0, n_tr), 
            rep(-sqrt(n_tr) + 1, n_tr))
  # computing beta
  cat("computing beta... \n")
  beta = solve.QP(Ktrtr, kappa, Amat, bvec)
  return(beta$solution)
}

### example
w <- kmm(as.matrix(sG_eur), as.matrix(sG_afr))
kmm_logreg <- cv.glmnet(sG_eur, sP_eur, family = "binomial",
                        weights = w) # error occured
