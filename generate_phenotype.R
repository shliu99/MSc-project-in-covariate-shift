# This file contains two functions for generating the phenotype for the
# European and African data 

# function 1: sample \beta
# @param: 
# len - total number of \beta being sampled
# alpha - proportion of causal features 
# seed - set seeds
# return: 
# beta - a matrix of one column which contains all the beta
# causal_ind - index of causal features

sample_beta <- function(len, alpha, seed = 2000){
  causal = alpha * len
  set.seed(seed)
  causal_ind = sample(1:len, size = causal)
  out = numeric(len)
  set.seed(seed + 20)
  out[causal_ind] = rnorm(causal)
  return(list(beta = as.matrix(out, ncol = 1), causal_ind = causal_ind))
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

generate_phenotype <- function(eur,afr,alpha,prop,upper,lower,seed = 1999){
  rate = TRUE
  count = 1
  n_eur = nrow(eur)
  n_afr = nrow(afr)
  p = ifelse(ncol(G_eur) == ncol(G_afr), ncol(G_eur), 
              stop("feature sizes do not match"))
  
  while (rate) {
    # generate beta
    b = sample_beta(p, alpha, seed + 10*count)
    beta = b$beta
    # generate phenotype
    linear_preds_eur = as.numeric(eur %*% beta)
    linear_preds_afr = as.numeric(afr %*% beta)
    sd_noise = sqrt(var(linear_preds_eur))
    set.seed(seed + 20*count)
    eps = rnorm(n_eur, mean = 0, sd = sd_noise)
    P_eur = linear_preds_eur + eps
    # now for the African
    set.seed(seed + count)
    P_afr = linear_preds_afr + rnorm(n_afr, mean = 0, sd = sd_noise)
    q = quantile(c(P_eur,P_afr), probs = (1 - prop))
    P_eur <- ifelse(P_eur >= q, 1, 0)
    P_afr <- ifelse(P_afr >= q, 1, 0)
    P_afr_mean = mean(P_afr); P_eur_mean = mean(P_eur)
    cat(paste("current proprtions of positve penotypes: European - ",
              P_eur_mean,'; African - ',P_afr_mean, '...\n'))
    if (P_afr_mean <= upper & P_afr_mean >= lower & 
        P_eur_mean <= upper & P_eur_mean >= lower) {rate = FALSE}
    count = count + 1
  }
  return(list(beta = beta, causal_ind = b$causal_ind,
              P_eur = P_eur, P_afr = P_afr))
}

