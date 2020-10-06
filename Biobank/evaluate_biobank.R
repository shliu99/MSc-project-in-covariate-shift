## evaluate UK biobank results 
library(readr)
library(pROC)

sdr = "spca"
dre = "ulsif"
weights_algo = paste(sdr, dre, sep = "_")

results <- matrix(NA, ncol = 7)
results = as.data.frame(results)
colnames(results) = c("algorithm", "fold", "auc", "threshold",
                      "sensitivity", "specificity", "accuracy")

for (fold in 1:5){
  res.path = file.path(paste0("fold", fold),
                       paste0(weights_algo,"_screen0.05_res.tsv"))
  preds = read_tsv(res.path)
  auc <- auc(preds$y, preds$pred_bin)
  # partial_auc <- auc(preds$y, preds$pred_bin, partial.auc = c(0,0.1),
  #                   partial.auc.correct = TRUE)
  roc.info <- roc(preds$y, preds$pred_bin)
  roc.df <- data.frame(sensitivity= roc.info$sensitivities,
                       specifity = roc.info$specificities,
                       threshold = roc.info$thresholds)
  roc.df$youden<- roc.df$sensitivity + roc.df$specifity - 1
  threshold.vals = roc.df[which(roc.df$youden == max(roc.df$youden)),]
  threshold <- threshold.vals$threshold
  labs <- ifelse(preds$pred_bin >= threshold, TRUE, FALSE)
  accurary <- sum(labs == preds$y)/nrow(preds)
  
  results <- rbind(results, c(weights_algo,
                              fold,
                              round(c(auc,
                              threshold,
                              threshold.vals$sensitivity,
                              threshold.vals$specifity,
                              accurary),3)))
}

results = results[-1,]
weights_algo
lapply(results[,c(3:7)], function(x) {round(mean(as.numeric(unlist(x))),3)})
lapply(results[,c(3:7)], function(x) {round(sd(as.numeric(unlist(x))),3)})
