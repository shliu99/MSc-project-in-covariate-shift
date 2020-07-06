# function that shows the performance of the prediction algorithm
# @param
# name: name of the functiion
# algo: a glmnet object
# test_data: design matrix of the test set data
# true: true labels 
# lambda: regularisation parameter
# threshold: threshold for converting the continuous predictions to binary ones 
# return:
# binary predictions
# coefficients of the algorithm


evaluate_algorithm <- function(name = "algorithm",algo,test_data,true,
                               lambda = "lambda.1se",
                               threshold = 0.5){
  n <- length(true)
  preds <- ifelse(predict(algo, test_data, s = lambda) >= threshold, 
                  1, 0)
  TP <-  sum(preds == 1 & true == 1)
  FP <- sum(preds == 1 & true == 0)
  TN <- sum(preds == 0 & true == 0)
  FN <- sum(preds == 0 & true == 1)
  # evaluation metrics
  cat(paste("\n",name,"\n"))
  cat(paste("\nNumber of parameters: ",sum(coef(algo) != 0)))
  cat(paste("\nPositive preds: "), sum(preds), "; True positives: ",sum(true))
  cat(paste("\nAccuracy: "), (TP+TN)/n)
  cat(paste("\nTrue positive rate: "), TP/(TP+FN))
  cat(paste("\nTrue negative rate: "), TN/(TN + FP))
  cat(paste("\nPrecision: ", TP/(TP + FP)),"\n")
  # roc
  cvec <- seq(0.001, 0.999, length.out = 1000)
  len <- length(cvec)
  TPR <- numeric(len)
  FPR <- numeric(len)
  for (k in 1:len) {
    TPR[k] <- sum(preds >= cvec[k] & true == 1)/sum(true == 1)
    FPR[k] <- sum(preds >= cvec[k] & true == 0)/sum(true == 0)
  }
  plot(y = TPR, x= FPR, type = "l", lwd = 2, col = "red",
       xlim = c(0,1), ylim = c(0,1), main = "ROC curve",
       ylab = "True positive rate", xlab = "False positive rate")
  lines(seq(0.001,0.999,length.out = 1000), seq(0.001,0.999,length.out = 1000))
  return(list(preds = preds, coeff = coef(algo), TPR = TPR, FPR = FPR))
}
