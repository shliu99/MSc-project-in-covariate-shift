library("mvtnorm"); library("ks"); library("AUC")
# X_tr = 0.75 N([0,...,0]^{d}, sigma) + 0.25 N([0.5,...,0.5]^{d}, sigma)
# x_te = 0.5 N([1.5,...,1.5]^{d}, sigma) + 0.5 N([2,...,2]^{d}, sigma)
# P(Y = 1 | X) = sigmoid(3X1 + sin(X2))

n = 250

lst <- list()

count = 1
for(count in 1:50) {
  results <- matrix(nrow = 19, ncol = 9)
  colnames(results) <- c(
    "kmml2",
    "kliepl2",
    "ulsifl2",
    "kdel2",
    "kmm%",
    "kliep%",
    "ulsif%",
    "kde%",
    "unweighted%"
  )
  rownames(results) <- paste0("d", 2:20)
  
  for (d in 2:20) {
    X.tr = matrix(NA, nrow = n, ncol = d)
    
    for (i in 1:n) {
      rand.num <- runif(1)
      if (rand.num <= 0.75) {
        X.tr[i, ] <- rmvnorm(n = 1,
                             mean = rep(0, d),
                             sigma = diag(d))
      } else{
        X.tr[i, ] <- rmvnorm(n = 1,
                             mean = rep(0.5, d),
                             sigma = diag(d))
      }
    }
    
    X.tr.pdf <-
      0.75 * dmvnorm(X.tr, mean = rep(0, d), sigma = diag(d)) +
      0.25 * dmvnorm(X.tr, mean = rep(0.5, d), sigma = diag(d))
    
    y.tr <- 1 / (1 + exp(-(3 * X.tr[, 1] + sin(X.tr[, 2]))))
    ##
    X.te = matrix(NA, nrow = n, ncol = d)
    
    for (i in 1:n) {
      rand.num <- runif(1)
      if (rand.num <= 0.5) {
        X.te[i, ] <- rmvnorm(n = 1,
                             mean = rep(1.5, d),
                             sigma = diag(d))
      } else{
        X.te[i, ] <- rmvnorm(n = 1,
                             mean = rep(2, d),
                             sigma = diag(d))
      }
    }
    
    X.te.pdf <-
      0.5 * dmvnorm(X.te, mean = rep(1.5, d), sigma = diag(d)) +
      0.5 * dmvnorm(X.te, mean = rep(2, d), sigma = diag(d))
    
    y.te <- 1 / (1 + exp(-(3 * X.te[, 1] + sin(X.te[, 2]))))
    
    ##
    true.ratio = X.tr.pdf / X.te.pdf
    kmm.ratio = kmm(X.tr, X.te)
    kliep.ratio = densratio(X.te, X.tr, method = "KLIEP")$compute_density_ratio(X.tr)
    ulsif.ratio = densratio(X.te, X.tr, method = "uLSIF")$compute_density_ratio(X.tr)
    if (d <= 0) {
      kde.ratio = predict(kde(X.tr), x = X.tr) / predict(kde(X.te), x = X.te)
    }
    
    results[d - 1, 1] <- sum((true.ratio - kmm.ratio) ^ 2)
    results[d - 1, 2] <- sum((true.ratio - kliep.ratio) ^ 2)
    results[d - 1, 3] <- sum((true.ratio - ulsif.ratio) ^ 2)
    if (d <= 0) {
      results[d - 1, 4] <- sum((true.ratio - kde.ratio) ^ 2)
    }
    
    y.tr <- ifelse(y.tr >= 0.5, 1, 0)
    y.te <- ifelse(y.te >= 0.5, 1, 0)
    
    # cv = cv.glmnet(x = X.tr, y = y.tr, family = "binomial", weights = kmm.ratio)
    # moda = glmnet(x = X.tr, y = y.tr, family = "binomial", weights = kmm.ratio,
    #               lambda = cv$lambda.1se)
    # predsa <- predict(moda, X.te)
    # results[, 5] <- sum((predsa - y.te)^2)
    #
    # cv = cv.glmnet(x = X.tr, y = y.tr, family = "binomial", weights = kliep.ratio)
    # modb = glmnet(x = X.tr, y = y.tr, family = "binomial", weights = kliep.ratio,
    #               lambda = cv$lambda.1se)
    # predsb <- predict(modb, X.te)
    # results[, 6] <- sum((predsb - y.te)^2)
    #
    # cv = cv.glmnet(x = X.tr, y = y.tr, family = "binomial", weights = ulsif.ratio)
    # modc = glmnet(x = X.tr, y = y.tr, family = "binomial", weights = ulsif.ratio,
    #               lambda = cv$lambda.1se)
    # predsc <- predict(modc, X.te)
    # results[, 7] <- sum((predsc - y.te)^2)
    #
    # if(d <= 4){
    #    cv = cv.glmnet(x = X.tr, y = y.tr, family = "binomial", weights = kde.ratio)
    #    modd = glmnet(x = X.tr, y = y.tr, family = "binomial", weights = kde.ratio,
    #                  lambda = cv$lambda.1se)
    #    predsd <- predict(modd, X.te)
    #    results[, 8] <- sum((predsd - y.te)^2)
    # }
    
    a = glm(y.tr ~ X.tr, family = "binomial", weights = kmm.ratio)
    predsa = predict(a, newdata = as.data.frame(X.te), type = "response")
    aa = roc(predsa, as.factor(y.te))
    results[d - 1, 5] <- auc(aa)
    
    b = glm(y.tr ~ X.tr, family = "binomial", weights = kliep.ratio)
    predsb = predict(b, newdata = as.data.frame(X.te), type = "response")
    bb = roc(predsb, as.factor(y.te))
    results[d - 1, 6] <- auc(bb)
    
    c = glm(y.tr ~ X.tr, family = "binomial", weights = ulsif.ratio)
    predsc = predict(c, newdata = as.data.frame(X.te), type = "response")
    cc = roc(predsc, as.factor(y.te))
    results[d - 1, 7] <- auc(cc)
    
    if (d <= 0) {
      e = glm(y.tr ~ X.tr, family = "binomial", weights = kde.ratio)
      predse = predict(e, newdata = as.data.frame(X.te), type = "response")
      ee = roc(predse, as.factor(y.te))
      results[d - 1, 8] <- auc(ee)
    }
    
    f = glm(y.tr ~ X.tr, family = "binomial")
    predsf = predict(f, newdata = as.data.frame(X.te), type = "response")
    ff = roc(predsf, as.factor(y.te))
    results[d - 1, 9] <- auc(ff)
  }
  lst[[count]] = results
  count = count + 1
}

produce_summary <- function(row, col, lst, out_mean, out_sd){
  temp <- vector("numeric", 50)
  for(i in 1:50){
    temp[i] = lst[[i]][row, col]
  }
  out_mean[row, col] <- mean(temp)
  out_sd[row, col] <- sd(temp)
  return(list(mean = out_mean, sd = out_sd))
}

out_mean <- matrix(NA, nrow = 19, ncol = 3)
rownames(out_mean) <- c("kmm","kliep","ulsif")
out_sd <- matrix(NA, nrow = 19, ncol = 3)
rownames(out_sd) <- c("kmm","kliep","ulsif")

for (row in 1:19) {
  for(col in 1:3){
    temp = vector("numeric")
    res = produce_summary(row, col, lst = lst, out_mean, out_sd)
    
  }
}

mean_results = matrix(nrow = 19, ncol = 9)
colnames(mean_results) <- c(
  "kmml2",
  "kliepl2",
  "ulsifl2",
  "kdel2",
  "kmm%",
  "kliep%",
  "ulsif%",
  "kde%",
  "unweighted%"
)
rownames(mean_results) <- paste0("d", 2:20)

sd_results = matrix(nrow = 19, ncol = 9)
colnames(sd_results) <- c(
  "kmml2",
  "kliepl2",
  "ulsifl2",
  "kdel2",
  "kmm%",
  "kliep%",
  "ulsif%",
  "kde%",
  "unweighted%"
)
rownames(sd_results) <- paste0("d", 2:20)


for (col in c(1:3, 5:7,9)){
  for(row in 1:19){
    temp = vector("numeric")
    for (elem in 1:50){
      temp = append(temp, lst[[elem]][row, col])
    }
    mean_results[row, col] <- mean(temp)
    sd_results[row, col] <- sd(temp)
  }
}


plot(x = 2:20, y = log10(mean_results[,1]), 
     xlab = "dimension", ylab = "L2 error of density ratio estimates (log10)",
     lwd = 2, type = "l", ylim = c(0,20))
lines(x = 2:20, y = log10(mean_results[,2]), col = "#FDAE61", lwd = 2)
lines(x = 2:20, y = log10(mean_results[,3]), col = "#4575B4", lwd = 2)
legend("bottomright", lty = c("solid", "solid", "solid"), 
       legend = c("KMM","KLIEP","uLSIF"),
       lwd = 2, col = c("black","#FDAE61", "#4575B4"), 
       bty = "n")


plot(x = 2:20, y = (mean_results[,5]), 
     xlab = "dimension", ylab = "auc",
     lwd = 2, type = "l")
lines(x = 2:20, y = (mean_results[,6]), col = "#FDAE61", lwd = 2)
lines(x = 2:20, y = (mean_results[,7]), col = "#4575B4", lwd = 2)
lines(x = 2:20, y = (mean_results[,9]), col = "#F46D43", lwd = 2)
legend("topright", lty = c("solid", "solid", "solid", "solid", "solid"), 
       legend = c("KMM","KLIEP","uLSIF", "unweighted"),
       lwd = 2, col = c("black","#FDAE61", "#4575B4", "#F46D43"), 
       bty = "n")
