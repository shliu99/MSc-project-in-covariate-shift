# install.packages("Matrix")
# install.packages("glmnet", repos = "http://cran.us.r-project.org")
# install.packages("densratio")
library("densratio"); library(glmnet); library("Matrix")
source("generating_phenotype.R")
source("subdata.R")
source("evaluate_algorithm.R")

# laod data ---------------------------------------------------------------
G_eur <- readMM("G_eur.mtx")
G_afr <- readMM("G_afr.mtx")
n_eur = nrow(G_eur)
n_afr = nrow(G_afr)
p <- ifelse(ncol(G_eur) == ncol(G_afr), ncol(G_eur), 
       stop("feature sizes do not match"))

# generate phenotype ------------------------------------------------------
# proprtion of positve African penotype is around 0.68
x <- generate_phenotype(G_eur, G_afr, alpha = 0.2, prop = 0.2, upper = 0.25, 
                        lower = 0.1)
beta <- x$beta; causal_ind <- x$causal_ind
P_afr <- x$P_afr; P_eur <- x$P_eur
rm(x)


# subset: sample size = 100, feature size = 200 ---------------------------
sub <- subdata(G_eur, P_eur, nsize = 100, psize = 200, causal = causal_ind)
sG_eur <- sub$data_sub; sP_eur <- sub$label_sub; scausal_ind <- sub$causal_sub
sub <- subdata(G_afr, P_afr, nsize = 100, psize = 200, causal = causal_ind)
sG_afr <- sub$data_sub; sP_afr <- sub$label_sub
rm(sub)
cat(paste("Number of postive European phenotype = ", sum(sP_eur),
    ";\nNumber of postive African phenotype = "), sum(sP_afr))

##########################################################################
# L1 binomial GLM
l1_logreg <- cv.glmnet(sG_eur, sP_eur, family = "binomial")
l1_logreg_res <- evaluate_algorithm("Logistic regression with L1 norm",
                                    l1_logreg, sG_afr, sP_afr)


