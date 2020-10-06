library("Matrix")

default_par <- par()
set.seed(1)
# load data ---------------------------------------------------------------
G_eur <- as.matrix(readMM("toy_data/G_eur.mtx"))
G_afr <- as.matrix(readMM("toy_data/G_afr.mtx"))

# simulate phenotypes -----------------------------------------------------
sample_beta <- function(len, alpha){
  causal = alpha * len
  causal_ind = sample(1:len, size = causal)
  out = numeric(len)
  out[causal_ind] = rnorm(causal)
  return(list(beta = as.matrix(out, ncol = 1), causal_ind = causal_ind))
}

generate_phenotype <- function(eur,afr,alpha,prop,upper,lower,seed = 1){
  rate = TRUE
  count = 1
  n_eur = nrow(eur)
  n_afr = nrow(afr)
  p = ifelse(ncol(G_eur) == ncol(G_afr), ncol(G_eur), 
             stop("feature sizes do not match"))
  
  while (rate) {
    # generate beta
    b = sample_beta(p, alpha)
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
    cat(paste("the current threshold for African samples is ", q, "\n"))
    if (P_afr_mean <= upper & P_afr_mean >= lower & 
        P_eur_mean <= upper & P_eur_mean >= lower) {rate = FALSE}
    count = count + 1
  }
  return(list(beta = beta, causal_ind = b$causal_ind,
              P_eur = P_eur, P_afr = P_afr))
}

x <- generate_phenotype(G_eur, G_afr, alpha = 0.2, prop = 0.2, upper = 0.25, 
                        lower = 0.1)
beta <- x$beta; causal_ind <- x$causal_ind
P_afr <- x$P_afr; P_eur <- x$P_eur

rm(x)

# write phenotypes 
write(P_eur, "eur_phe.txt", ncol = 1)
write(P_afr, "afr_phe.txt", ncol = 1)
