#!/apps/eb/skylake/software/R/3.6.2-foss-2019b/bin/Rscript
#$ -cwd
#$ -q short.qe
#$ -P mcvean.prjc
#$ -pe shmem 8
#$ -M wolf5659@ox.ac.uk
#$ -o logs/toy_weights.o.log
#$ -e logs/toy_weights.e.log
#$ -m ea


# compute weights for toy data

prop = 0.2 #ratio of black to white in each train fold
sdr = "kspca"

library("lfda"); library("SpNMF"); library("NMF");library("densratio")
library("quadprog");library("Matrix")


## prepare data set --------------------------------------------------------
data.X.train <- as.matrix(readMM("toy_data/G_eur.mtx"))
data.X.test <- as.matrix(readMM("toy_data/G_afr.mtx"))
data.y = as.numeric(readLines("eur_phe.txt"))
data.y.true = as.numeric(readLines("afr_phe.txt"))


## supervised dimension reduction ------------------------------------------

# kernel functions - rbf and linear kernels
source("scripts/fun_kernel.R")
# Gaussian width selection
source("scripts/fun_gaussian_cv.R")

# (1) PCA
# SPCA
source("scripts/fun_spca_dual.R")
# KSPCA
source("scripts/f_fun_kspca.R")

# (2) FDA
# lfda
source("scripts/fun_lfda.R")
# klfda
source("scripts/fun_klfda.R")

## density ratio estimation
# (1) KMM
source("scripts/fun_kmm.R")
# (2) supporting functions for KLIEP
source("scripts/fun_kliep.R")


gaussian_width <- kernelWidth_median(data.X.train)

sdr_lst <- c("lfda","klfda")  # "kspca", "lfda", "klfda" "snmf"
for (sdr in sdr_lst){
  switch (sdr,
          spca = {
            spca.res <- spca(data.X.train, data.y)
            X.train.out <- data.X.train %*% spca.res
            X.test.out <- data.X.test %*% spca.res
          },
          kspca = {
            kspca.res <- kSPCA(data.X.train, data.y, gaussian_width)
            X.train.out <- rbfkernel(data.X.train, sigma = gaussian_width) %*% kspca.res
            X.test.out <- rbfkernel(data.X.test, data.X.train, sigma = gaussian_width) %*% kspca.res
          },
          lfda = {
            lfda.res <- lfda(data.X.train, data.y, 1,
                             knn = 7 ,metric = "orthonormalized")
            X.train.out <- predict(lfda.res, data.X.train)
            X.test.out <- predict(lfda.res, data.X.test)
          },
          klfda = {
            klfda.res <- klfda(rbfkernel(data.X.train, sigma = gaussian_width), 
                               data.y, 
                               1, 
                               knn = 7 ,metric = "orthonormalized")
            X.train.out <- klfda.res$Z
            X.test.out <- rbfkernel(data.X.test, data.X.train, 
                                    sigma = gaussian_width) %*% klfda.res$T  
          },
          snmf = {
            data.X.train <- data.X.train[,colSums(data.X.train)!=0]
            data.X.test <- data.X.test[,colSums(data.X.test)!=0]
            t_types <-  chty(data.X.train, data.y, 2, 3)
            nmb1 <-  t_types$r1
            nmb2 <-  t_types$r2
            t_total <- getT(data.X.train,data.y,nmb1,nmb2)
            X.train.out <- spnmf(data.X.train, t_total)$W
            X.test.out <- spnmf(data.X.test, t_total)$W
          },
          default = {
            X.train.out <-  data.X.train
            X.test.out <- data.X.test
          }
  )
  write.csv(X.train.out, file = file.path(fold_path, paste0(sdr,"_train.csv")),
            row.names = F)
  write.csv(X.test.out, file = file.path(fold_path, paste0(sdr,"_test.csv")),
            row.names = F)
}




