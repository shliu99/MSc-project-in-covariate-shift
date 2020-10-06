#!/apps/eb/skylake/software/R/3.6.2-foss-2019b/bin/Rscript
#$ -cwd
#$ -P mcvean.prjc
#$ -N prescreen
#$ -pe shmem 8
#$ -q short.qe
#$ -o logs/prescreen.o.log
#$ -e logs/prescreen.e.log


# compute p values for each variant in each fold

plink2_path  <- "/apps/eb/skylake/software/PLINK/2.00a2.3_x86_64/plink2"
zstdcat_path <- "/apps/eb/skylake/software/zstd/1.4.4-GCCcore-8.3.0/bin/zstdcat"

# devtools::install_github("brieuclehmann/snpnet")
library(devtools)
library(snpnet)
library(readr)
library(dplyr)
library(parallel)
# install_github("junyangq/glmnetPlus")
# install_github("chrchang/plink-ng", subdir="/2.0/cindex")
# install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
set.seed(1)

pheno  <- "breast"
prop   <- 0.2
nfolds <- 5
cl = parallel::makeCluster(detectCores())

for (fold in 1:nfolds){
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
  # vars <- dplyr::mutate(dplyr::rename(
  #   data.table::fread(cmd=paste0(configs[['zstdcat.path']], ' ',
  #                                paste0(pfile, '.pvar.zst'))), 'CHROM'='#CHROM'),
  #   VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  # n_genotypes <- length(vars)
  
  # for an R interactive session
  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(zstdcat_path, ' ', paste0(gfile, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  n_genotypes <- length(vars)
  
  
  
  pvar <- pgenlibr::NewPvar(paste0(gfile, '.pvar.zst'))
  pgen <- pgenlibr::NewPgen(paste0(gfile, '.pgen'), pvar = pvar,
                            sample_subset = match(phe$ID, ids))
  data.y <- phe[[pheno]]
  # genotype data
  if(pheno == "breast" | pheno == "prostate"){
    data.X.gen2 <- as.matrix(cbind(data.y, Age = phe$Age,
                                   phe[, paste0("PC", 1:10)]))
  } else{
    data.X.gen2 <- as.matrix(cbind(data.y, Age = phe$Age, Sex1 = ifelse(phe$Sex==1,1,0),
                                   Sex0 = ifelse(phe$Sex == 0,1,0), 
                                   phe[, paste0("PC", 1:10)], phe[, paste0("PC", 1:10, "_Sex")]))
  }
  
  
  # fit lm on variant ~ phe + sex + age + pc + pc*sex
  p = ncol(data.X.gen2)
  n = nrow(data.X.gen2)
  XtX.inv = solve(t(data.X.gen2) %*% data.X.gen2)
  XtX.inv.Xt =  XtX.inv %*% t(data.X.gen2)
  
  
  prescreening <- function(col, data.X.gen2, p, n, XtX.inv, XtX.inv.Xt){
    col = matrix(col, ncol = 1)
    beta.hat = matrix(XtX.inv.Xt %*% col, ncol = 1)
    RSS = t(col - data.X.gen2%*%beta.hat) %*% (col - data.X.gen2%*%beta.hat)
    sigma2.hat = as.numeric(RSS / (n - p))
    beta.var = XtX.inv[1,1]*sigma2.hat
    p.val = 2*pt(abs(beta.hat[1]/sqrt(beta.var)), df = (n-p), lower.tail = F)
    return(p.val)
  }
  
  p.vals = vector("numeric")
  
  ### parallel computing
  count = 0; n_each_time = 5000
  
  while (count +  n_each_time <= n_genotypes){
    data.X.train <- pgenlibr::ReadList(pgen, 1:5000, meanimpute=T)
    p.vals = append(p.vals, parallel::parApply(cl,
                                               data.X.train, MARGIN = 2, prescreening, 
                                               data.X.gen2 = data.X.gen2, p = p, n = n,
                                               XtX.inv = XtX.inv, XtX.inv.Xt = XtX.inv.Xt))
    count <- count +  n_each_time
    cat(paste("current sample: ", count))
  }
  
  if(count +  n_each_time > n_genotypes) {
    data.X.train <- pgenlibr::ReadList(pgen, (count+1):n_genotypes, meanimpute=T) 
    p.vals = append(p.vals, parallel::parApply(cl,
                                               data.X.train, MARGIN = 2, prescreening, 
                                               data.X.gen2 = data.X.gen2, p = p, n = n,
                                               XtX.inv = XtX.inv, XtX.inv.Xt = XtX.inv.Xt))}
  
  # parallel::stopCluster(cl)
  
  features_selected = matrix(NA, ncol = 2, nrow = n_genotypes)
  colnames(features_selected) = c("inx","p.val")
  features_selected[,1] = 1:n_genotypes
  features_selected[,2] = p.vals
  features_selected <- features_selected[order(features_selected[, 2]),]
  outfile = file.path("data", pheno,
                      paste0("fold", fold),
                      paste0("prop", prop, "_pvals.txt"))
  write(t(features_selected), outfile, ncolumns  = 2)
  rm(p.vals,features_selected)
}