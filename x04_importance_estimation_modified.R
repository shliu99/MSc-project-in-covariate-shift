#!/apps/eb/skylake/software/R/3.6.2-foss-2019b/bin/Rscript
#$ -cwd
#$ -q short.qe
#$ -P mcvean.prjc
#$ -pe shmem 8
#$ -o logs/x04_estimate_weights.o.log
#$ -e logs/x04_estimate_weights.e.log
# This estimate weights of the training samples

### Set up parameters ---

plink2_path  <- "/apps/eb/skylake/software/PLINK/2.00a2.3_x86_64/plink2"
zstdcat_path <- "/apps/eb/skylake/software/zstd/1.4.4-GCCcore-8.3.0/bin/zstdcat"

# devtools::install_github("brieuclehmann/snpnet")
library(devtools)
library(snpnet)
library(readr)
library(dplyr)

library("lfda"); library("SpNMF"); library("NMF")
library("densratio"); library("quadprog"); library("Matrix")

# install_github("junyangq/glmnetPlus")
# install_github("chrchang/plink-ng", subdir="/2.0/cindex")
# install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
set.seed(1)

ncores <- as.integer(Sys.getenv("NSLOTS"))
pheno <- Sys.getenv("pheno")
prop <- as.double(Sys.getenv("prop_black"))
prop_keep <- as.double(Sys.getenv("prop_screen"))
fold <- as.integer(Sys.getenv("fold"))
sdr <- Sys.getenv("sdr")
dre <- Sys.getenv("dre")
nfolds <- 5
weights_algo = paste(sdr, dre, sep = "_")
# 
# pheno  <- "breast"
# prop   <- 0.2
# fold   <- 1
# nfolds <- 5
# prop_keep = 0.05

kfile  <- file.path("data", pheno,
                    paste0("fold", fold),
                    paste0("prop", prop, ".txt"))
file_test <- file.path("data", pheno,
                       paste0("fold", fold),
                       paste0("prop", prop, "_test.txt"))
pfile  <- "data/all_vars.tsv"
gfile  <- "data/ukb_cal_pgen/ukb_cal_all"
covars <- c("Sex", "Age", paste0("PC", 1:10), paste0("PC", 1:10, "_Sex"))
configs <- list(results.dir = file.path("temp", pheno,
                                        paste0("fold", fold),
                                        paste0("prop", prop)),
                plink2.path = plink2_path,
                zstdcat.path = zstdcat_path,
                keep = kfile)

train_ids <- read.table(kfile)[, 1]
ntrain <- length(train_ids)
test_ids <- read.table(file_test)[,1]
ntest <- length(test_ids)

ids <- readIDsFromPsam(paste0(gfile, ".psam"))
id_df <- data.frame(ID = ids, stringsAsFactors = F) %>%
  mutate(sort_order = 1:n())

### pre-screening

phe <- readPheMaster(phenotype.file = pfile, psam.ids = ids,
                     family = NULL, covariates = covars,
                     phenotype = pheno, status = NULL,
                     split.col = NULL, configs = configs) # in configs, keep = train_fold file
# for shell command
vars <- dplyr::mutate(dplyr::rename(
  data.table::fread(cmd=paste0(configs[['zstdcat.path']], ' ',
                               paste0(gfile, '.pvar.zst'))), 'CHROM'='#CHROM'),
  VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID

# for an R interactive session
# vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(zstdcat_path, ' ', paste0(gfile, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
# n_genotypes <- length(vars)

n_genotypes <- length(vars)


pvar <- pgenlibr::NewPvar(paste0(gfile, '.pvar.zst'))
pgen <- pgenlibr::NewPgen(paste0(gfile, '.pgen'), pvar = pvar,
                          sample_subset = match(phe$ID, ids))
data.y <- phe[[pheno]]

# genotype data - train
p.vals.path = file.path("data", pheno,
                        paste0("fold", fold),
                        paste0("prop", prop, "_pvals.txt"))
p.vals = read.table(p.vals.path, header = F)
inx = p.vals[1:floor(prop_keep*n_genotypes),1]

# genotype data - train
# data.X.train = pgenlibr::ReadList(pgen, inx, meanimpute=T)
inx = inx[order(inx)]
spt = split(seq_len(n_genotypes),cut(seq_len(n_genotypes),20))

data.X.temp = pgenlibr::ReadList(pgen, spt[[1]], meanimpute=F)
inx.temp = which(spt[[1]] %in% inx)
data.X.train = data.X.temp[,inx.temp]
for(i in 2:20){
  data.X.temp = pgenlibr::ReadList(pgen, spt[[i]], meanimpute=F)
  inx.temp = which(spt[[i]] %in% inx)
  data.X.train = cbind(data.X.train, data.X.temp[,inx.temp])
}
rm(data.X.temp, inx.temp)
colnames(data.X.train) = vars[inx]


# genotype data - test
# data.X.test <- pgenlibr::ReadList(pgen, inx, meanimpute=T)
configs <- list(results.dir = file.path("temp", pheno,
                                        paste0("fold", fold),
                                        paste0("prop", prop)),
                plink2.path = plink2_path,
                zstdcat.path = zstdcat_path,
                keep = file_test)
phe <- readPheMaster(phenotype.file = pfile, psam.ids = ids,
                     family = NULL, covariates = covars,
                     phenotype = pheno, status = NULL,
                     split.col = NULL, configs = configs) # configs keep = test fold
pgen <- pgenlibr::NewPgen(paste0(gfile, '.pgen'), pvar = pvar,
                          sample_subset = match(phe$ID, ids))

data.X.temp = pgenlibr::ReadList(pgen, spt[[1]], meanimpute=F)
inx.temp = which(spt[[1]] %in% inx)
data.X.test = data.X.temp[,inx.temp]
for(i in 2:20){
  data.X.temp = pgenlibr::ReadList(pgen, spt[[i]], meanimpute=F)
  inx.temp = which(spt[[i]] %in% inx)
  data.X.test = cbind(data.X.test, data.X.temp[,inx.temp])
}
rm(data.X.temp, inx.temp)
colnames(data.X.test) = vars[inx]

# mean imputation and remove columns that are all NAs
p <- length(inx)
pnas_train <- numeric(p)
pnas_test <- numeric(p)
for (j in 1:p) {
  pnas_train[j] <- mean(is.na(data.X.train[, j]))
  data.X.train[is.na(data.X.train[, j]), j] <- mean(data.X.train[, j], na.rm = T) 
  pnas_test[j] <- mean(is.na(data.X.test[, j]))
  data.X.test[is.na(data.X.test[, j]), j] <- mean(data.X.test[, j], na.rm = T) 
}

all.nas = unique(c(which(pnas_train==1), which(pnas_test==1)))

# remove rows that are duplicated and remove columns that are all NAs
data.X.train = unique(data.X.train[,!(1:p %in% all.nas)], MARGIN = 1)
data.X.test = unique(data.X.test[,!(1:p %in% all.nas)], MARGIN = 1)

rm(phe)
# now we have data.X.train, data.X.test, data.y ready




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


## Supervised dimension reduction - implementation
switch (pheno,
  breast = {gaussian_width = 7400},
  prostate = {gaussian_width = 550}
) # (from 03d cv_on_rbf)


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

rm(data.X.train, data.X.test)

cat(paste("the current phenotype is ", pheno, " prop = ", prop))
cat(paste("the current algorithm is: ",weights_algo,"\n"))
cat(paste("For fold ",fold," with ", sdr,
          "the dimension of reduced space is ",ncol(X.train.out),"\n"))

## density ratio estimation - implementation
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
outfile <- file.path("data", pheno,
                     paste0("fold", fold),
                     paste0("prop",prop,"_weights"))

for (f in seq_len(nfolds)) {
  dir.create(outfile,
             showWarnings = FALSE, recursive = TRUE)
}

write(weights, file = file.path(outfile,paste0(weights_algo,
                                               "_screen",prop_keep,
                                               ".txt")), ncolumns = 1)


sdr.path <- file.path("data", pheno,
                      paste0("fold", fold),
                      paste0("prop",prop))

dir.create(sdr.path, showWarnings = F, recursive = T)

write.csv(X.train.out, file = file.path(sdr.path, paste0(sdr,"_train.csv")),
          row.names = F)
write.csv(X.test.out, file = file.path(sdr.path, paste0(sdr,"_test.csv")),
          row.names = F)
