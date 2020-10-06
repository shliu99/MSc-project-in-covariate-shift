#!/apps/eb/skylake/software/R/3.6.2-foss-2019b/bin/Rscript
#$ -cwd
#$ -P mcvean.prjc
#$ -N lfda_cv
#$ -pe shmem 8
#$ -q short.qe
#$ -o logs/lfda_cv.o.log
#$ -e logs/lfda_cv.e.log

# experiment in choosing the dimensionality of lfda/klfda reduced space


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
pheno  <- "breast"
prop   <- 0.2
fold   <- 1
nfolds <- 5
prop_keep <- 0.05

kfile  <- file.path("data", pheno,
                    paste0("fold", fold),
                    paste0("prop", prop, ".txt"))
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
n_genotypes <- length(vars)

# for an R interactive session
# vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(zstdcat_path, ' ', paste0(gfile, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
# n_genotypes <- length(vars)

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
# colnames(data.X.train) = vars[inx]


# mean imputation and remove columns that are all NAs
p <- length(inx)
pnas_train <- numeric(p)
for (j in 1:p) {
  pnas_train[j] <- mean(is.na(data.X.train[, j]))
  data.X.train[is.na(data.X.train[, j]), j] <- mean(data.X.train[, j], na.rm = T) 
}

all.nas = which(pnas_train==1)
data.X.train = data.X.train[,!(1:p %in% all.nas)]

## cv for choosing the dimensionality
library("class");library("lfda")
source("scripts/fun_lfda.R")

N = 1000
set.seed(10)
data.cv = sample(1:ntrain, 1000)
data.X.train = data.X.train[data.cv,]
train.inx = list(); test.inx = list()
inx = split(seq_len(N), cut(seq_len(N), breaks = 5))
for(i in 1:5){
  train.inx[[i]] = seq_len(N)[!(seq_len(N) %in% inx[[i]])]
  test.inx[[i]] = inx[[i]]
}

dim.lst = c(10,100,300,500)
mis.class.error = matrix(NA, ncol = 5, nrow = length(dim.lst))
rownames(mis.class.error) = dim.lst
colnames(mis.class.error) = paste0("foldi",1:5)

for (foldi in 1:5) {
  data.X.train.cv = data.X.train[train.inx[[foldi]],]
  data.X.test.cv = data.X.train[test.inx[[foldi]],]
  data.y.cv = data.y[train.inx[[foldi]]]
  data.y.test = data.y[test.inx[[foldi]]]
  for(dim in dim.lst){
    lfda.res <- lfda(data.X.train.cv, data.y.cv, dim,
                     knn = 7 ,metric = "orthonormalized")
    X.train.lfda <- predict(lfda.res, data.X.train.cv)
    X.test.lfda <- predict(lfda.res, data.X.test.cv)
    knn.res = knn1(X.train.lfda, X.test.lfda, data.y.cv)
    mis.rate = sum(knn.res != data.y.test)/length(knn.res)
    mis.class.error[which(dim.lst %in% dim), foldi] = mis.rate
    cat(paste("the misclassification rate for cv fold ", foldi, "dimension ", dim, "is", mis.rate,"\n"))
  }
}


cat(paste("the dim.lst is: ", dim.lst,"\n"))
print(paste("the misclassification errors are: ", rowMeans(mis.class.error), "\n"))
write.table(mis.class.error, file.path(configs$results.dir,"_lfda.csv"))
