#!/apps/eb/skylake/software/R/3.6.2-foss-2019b/bin/Rscript
#$ -cwd
#$ -q short.qe
#$ -P mcvean.prjc
#$ -pe shmem 2
#$ -o logs/gaussian_width_cv.log
#$ -e logs/gaussian_width_cv.log
# This chooses the rbf kernel parameter


plink2_path  <- "/apps/eb/skylake/software/PLINK/2.00a2.3_x86_64/plink2"
zstdcat_path <- "/apps/eb/skylake/software/zstd/1.4.4-GCCcore-8.3.0/bin/zstdcat"

# devtools::install_github("brieuclehmann/snpnet")
library(devtools)
library(snpnet)
library(readr)
library(dplyr)
# install_github("junyangq/glmnetPlus")
# install_github("chrchang/plink-ng", subdir="/2.0/cindex")
# install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
set.seed(1)

ncores <- as.integer(Sys.getenv("NSLOTS"))
pheno  <- "prostate"
prop   <- 0.2
fold   <- 1
nfolds <- 5
prop_keep <- 0.05

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
# vars <- dplyr::mutate(dplyr::rename(
#   data.table::fread(cmd=paste0(configs[['zstdcat.path']], ' ',
#                                paste0(gfile, '.pvar.zst'))), 'CHROM'='#CHROM'),
#   VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
# n_genotypes <- length(vars)

# for an R interactive session
vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(zstdcat_path, ' ', paste0(gfile, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
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

p <- length(inx)
for (j in 1:p) {
  data.X.train[is.na(data.X.train[, j]), j] <- mean(data.X.train[, j], na.rm = T) 
}

### cv on gaussian width

rbfkernel_cv <- function(X, sigma_lst, Y = NULL, lab) {
  if (!is.matrix(X))
  {print("X must be a matrix containing samples in its rows")
    return()
  }
  if (!is.null(Y)){
    if (!is.matrix(Y))
    {print("Y must be a matrix containing samples in its rows or NULL if it should not be used")
      return()
    }
    if (ncol(X) != ncol(Y)){
      print("The samples in the rows of X and Y must be of same dimension")
      return()
    }
  }
  lab = as.numeric(lab) + 1
  n <- nrow(X)
  if (is.null(Y)) {
    XtX <- tcrossprod(X)
    XX <- matrix(1, n) %*% diag(XtX)
    D <- XX - 2 * XtX + t(XX)
  }
  else{
    m <- nrow(Y)
    XX <- matrix(apply(X ^ 2, 1, sum), n, m)
    YY <- matrix(apply(Y ^ 2, 1, sum), n, m, byrow = TRUE)
    XY <- tcrossprod(X, Y)
    D <- XX - 2 * XY + YY
  }
  # calculate rbf-kernel matrix
  cv.res = matrix(NA, nrow = 2, ncol = length(sigma_lst))
  cv.res[1,] = sigma_lst
  for (sigma in sigma_lst){
    K <-  exp(-D / (2 * sigma))
    withinSum <-  0
    num = matrix(0,1,max(lab))
    for (i in 0:max(lab)) {
      index <- lab == i
      withinSum <-  withinSum + sum(K[index, index])
      num[i] <-  sum(index)
    }
    betweenSum <-  sum(K) - withinSum
    withinNum <-  sum(num^2)
    betweenNum <-  sum(num)^2 - withinNum
    w <-  withinSum/withinNum
    b <-  betweenSum/betweenNum
    j <-  1 - w + b
    cv.res[2,which(sigma_lst %in% sigma)] <- j 
  }
  return(cv.res)
}


# res1 <- rbfkernel_cv(X = data.X.train, sigma_lst = c(1,5,10,20,50,100,1000,5000,10000,15000), lab = data.y)
# res2 <- rbfkernel_cv(X = data.X.train, sigma_lst = c(seq(100,1000,100)), lab = data.y)
# res3 <- rbfkernel_cv(X = data.X.train, sigma_lst = c(seq(500,600,10)), lab = data.y)
# res4 <- rbfkernel_cv(X = data.X.train, sigma_lst = c(seq(7000,8000,100)), lab = data.y)
# best sigma for "breast" ~= 7400; for "prostate" ~= 550

cat(paste("sigmas are: ", res[1,]))
cat("\n")
cat(paste("scores are: ", res[2,]))



