# install.packages("Matrix")
# install.packages("glmnet", repos = "http://cran.us.r-project.org")
# install.packages("densratio")
library("densratio"); library(glmnet); library("Matrix")
source("generating_phenotype.R")


# laod data ---------------------------------------------------------------
G_eur <- readMM("G_eur.mtx")
G_afr <- readMM("G_afr.mtx")
n_eur = nrow(G_eur)
n_afr = nrow(G_afr)
ifelse(ncol(G_eur) == ncol(G_afr), p <- ncol(G_eur), 
       stop("feature sizes do not match"))

# generate phenotype ------------------------------------------------------
x <- sample_beta(p, 0.05, seed = 2000)
beta <- x$beta; casual_ind <- x$casual_ind
y <- generate_phenotype(G_eur, G_afr, beta, prop = 0.2, upper = 0.25, 
                        lower = 0.1)
# proprtion of positve African penotype is around 0.68
x <- sample_beta(p, 0.05, seed = 20)
beta <- x$beta; casual_ind <- x$casual_ind
y <- generate_phenotype(G_eur, G_afr, beta, prop = 0.2, upper = 0.25, 
                        lower = 0.1)
# proprtion of positve African penotype is 0.22
x <- sample_beta(p, 0.05, seed = 100)
beta <- x$beta; casual_ind <- x$casual_ind
y <- generate_phenotype(G_eur, G_afr, beta, prop = 0.2, upper = 0.25, 
                        lower = 0.1)
# proprtion of positve African penotype is 0.16
x <- sample_beta(p, 0.05, seed = 999)
beta <- x$beta; casual_ind <- x$casual_ind
y <- generate_phenotype(G_eur, G_afr, beta, prop = 0.2, upper = 0.25, 
                        lower = 0.1)
# proprtion of positve African penotype is 0.18

G_afr <- y$G_afr; G_eur <- y$G_eur
