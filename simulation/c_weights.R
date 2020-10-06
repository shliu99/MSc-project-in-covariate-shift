## density ratio estimation - implementation
library(densratio)

sdr = "snmf"
dre = "kliep"
# weights_algo = paste(sdr, dre, sep = "_")

X.train.out <- read.csv(file = paste0(sdr,"_train.csv"))[,-1]
X.test.out <- read.csv(file = paste0(sdr,"_test.csv"))[,-1]

# kernel functions - rbf and linear kernels
source("scripts/fun_kernel.R")
# Gaussian width selection
source("scripts/fun_gaussian_cv.R")
## density ratio estimation
# (1) KMM
source("scripts/fun_kmm.R")
# (2) supporting functions for KLIEP
source("scripts/fun_kliep.R")

dre = "ulsif"

switch (dre,
        kmm = {
          weights <- kmm(X.train.out, X.test.out)
        },
        kliep = {
          num <- nrow(X.test.out) * 0.5
          kliep.obj<- densratio(X.test.out, X.train.out, method = "KLIEP", kernel_num = num)
          weights <- compute_kernel_Gaussian(X.train.out, kliep.obj$kernel_info$centers, kliep.obj$kernel_info$sigma) %*% kliep.obj$kernel_weights
        },
        ulsif = {
          num <- nrow(X.test.out) * 0.5
          weights <- densratio(X.test.out, X.train.out, method = "uLSIF", kernel_num = num)$compute_density_ratio(X.train.out)
        },
        default = {
          weights = rep(1, ntrain)
        }
)


## output results
weights_algo = paste(sdr, dre, sep = "_")
out_path = file.path("weights")
dir.create(out_path, showWarnings = F, recursive = T)
write(weights, file.path(out_path, paste0(weights_algo,".txt")), ncolumns = 1)
