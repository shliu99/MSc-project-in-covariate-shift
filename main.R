# install.packages("Matrix")
# install.packages("glmnet", repos = "http://cran.us.r-project.org")
# install.packages("densratio")
library("densratio"); library(glmnet); library("Matrix")
source("generating_phenotype.R")
source("subdata.R")

# laod data ---------------------------------------------------------------
G_eur <- readMM("G_eur.mtx")
G_afr <- readMM("G_afr.mtx")
n_eur = nrow(G_eur)
n_afr = nrow(G_afr)
ifelse(ncol(G_eur) == ncol(G_afr), p <- ncol(G_eur), 
       stop("feature sizes do not match"))

# generate phenotype ------------------------------------------------------
# proprtion of positve African penotype is around 0.68
x <- sample_beta(p, 0.05, seed = 20)
beta <- x$beta; causal_ind <- x$causal_ind
y <- generate_phenotype(G_eur, G_afr, beta, prop = 0.2, upper = 0.25, 
                        lower = 0.1)
P_afr <- y$P_afr; P_eur <- y$P_eur
rm(x,y)



