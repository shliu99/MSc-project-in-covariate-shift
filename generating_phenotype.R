# This file contains two functions for generating the phenotype for the
# European and African data 

# function 1: sample \beta
# @param: 
# len - total number of \beta being sampled
# alpha - proportion of casual features 
# seed - set seeds
# return: 
# beta - a matrix of one column which contains all the beta
# casual_ind - index of causual features

sample_beta <- function(len, alpha, seed = 2000){
  set.seed(seed)
  casual = alpha * len
  casual_ind = sample(1:len, size = casual)
  out = numeric(len)
  set.seed(seed + 20)
  out[casual_ind] = rnorm(casual)
  return(list(beta = as.matrix(out, ncol = 1), casual_ind = casual_ind))
}

# function 2: generate phenotye for the European and African
# @param
# eur: the design matrix containing the European data
# afr: the design matrix containing the Afircan data
# beta: beta simulated from function sample_beta
# prop: proportion of Eur phenotype assigned to +1
# lower: the minimal threshold for the porportion of positive phenotype in Afr
# upper: the maximal threshold for the porportion of positive phenotype in Afr
# seed - set seeds

generate_phenotype <- function(eur,afr,beta,prop,upper,lower,seed = 1999){
  rate = TRUE
  n_eur = nrow(eur)
  n_afr = nrow(afr)
  linear_preds_eur = as.numeric(eur %*% beta)
  linear_preds_afr = as.numeric(afr %*% beta)
  sd_noise = sqrt(var(linear_preds_eur))
  count = 1
    while (rate) {
    count = count + 10
    set.seed(seed + count)
    eps = rnorm(n_eur, mean = 0, sd = sd_noise)
    P_eur = linear_preds_eur + eps
    q = quantile(P_eur, probs = (1 - prop))
    P_eur <- ifelse(P_eur >= q, 1, 0)
    # now for the African
    set.seed(seed + count)
    P_afr = linear_preds_afr + rnorm(n_afr, mean = 0, sd = sd_noise)
    P_afr <- ifelse(P_afr >= q, 1, 0)
    P_afr_mean = mean(P_afr)
    cat(paste("current proprtion of positve African penotype is ",
              P_afr_mean,'...\n'))
    if (P_afr_mean <= upper & P_afr_mean >= lower) {rate = FALSE}
  }
  return(list(P_eur = P_eur, P_afr = P_afr))
}

