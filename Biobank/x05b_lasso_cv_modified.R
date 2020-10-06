#!/apps/eb/skylake/software/R/3.6.2-foss-2019b/bin/Rscript
#$ -cwd
#$ -q short.qe
#$ -P mcvean.prjc
#$ -pe shmem 6
#$ -o logs/x05b_lasso_cv.o.log
#$ -e logs/x05b_lasso_cv.e.log

### Set up parameters ---

# devtools::install_github("brieuclehmann/snpnet")
library(snpnet)
library(readr)
library(dplyr)
set.seed(1)

ncores <- as.integer(Sys.getenv("NSLOTS"))
mem <- ncores * 10e3
pheno <- Sys.getenv("pheno")
prop <- as.double(Sys.getenv("prop_black"))
prop_keep <- as.double(Sys.getenv("prop_screen"))
fold <- as.integer(Sys.getenv("fold"))
i <- as.integer(Sys.getenv("SGE_TASK_ID"))
nfolds <- 5
sdr <- Sys.getenv("sdr")
dre <- Sys.getenv("dre")
weights_algo = paste(sdr, dre, sep = "_")
gamma = 0.5

# 
# 
# ncores <- 24
# mem <- ncores * 10e3
# pheno  <- "breast"
# prop   <- 0.2
# fold   <- 1
# nfolds <- 5
# prop_keep = 0.05
# i <- 1
# sdr = "lfda" # spca kspca lfda klfda snmf
# dre = "ulsif" # kmm kliep ulsif
# weights_algo = paste(sdr, dre, sep = "_")
# gamma = 0.5

kfile  <- file.path("data", pheno,
                    paste0("fold", fold),
                    paste0("prop", prop, ".txt"))
pfile <- "data/all_vars.tsv"
gfile <- "data/ukb_cal_pgen/ukb_cal_all"

results.dir <- file.path("temp", pheno,
                         paste0("fold", fold),
                         paste0("prop", prop))
configs <- readRDS(file.path(results.dir, 
                             paste0(weights_algo, "_screen",prop_keep,"_flat_configs.RDS")))
full.lams <- scan(file.path(configs$results.dir, paste0(weights_algo,
                                                        "_screen",prop_keep,
                                                        "_flat_full.lams.txt")))
if (pheno %in% c("breast", "cancer")) {
  covars <- c("Age", paste0("PC", 1:10))
} else {
  covars <- c("Sex", "Age", paste0("PC", 1:10), paste0("PC", 1:10, "_Sex"))
}

cv.phenotype.file <- paste0(configs[["cv.full.prefix"]], ".tsv")
################################################################################
weights_path <- file.path("data", pheno,
                          paste0("fold", fold),
                          paste0("prop",prop,"_weights"))
weights <- as.numeric(readLines(file.path(weights_path,paste0("default_default",
                                                         "_screen",prop_keep,
                                                         ".txt"))))
weights <- weights^gamma
weights = weights/mean(weights)
weights[weights > 5] <- 5
weights[weights <= 1e-3] <- 0
################################################################################

cv_configs <- configs
cv_configs$results.dir <- paste0(configs$results.dir, "/foldi", i)
cv_configs$nCores <- ncores
cv_configs$mem <- mem
cv_configs$stopping.lag <- 5
cv_configs$gcount.full.prefix <- NULL


dir.create(cv_configs$results.dir,
             showWarnings = FALSE, recursive = TRUE)

### Fit lasso for ith fold ---

out <- snpnet(gfile,
              cv.phenotype.file,
              pheno,
              covariates = covars,
              weights = weights,
              full.lams = full.lams,
              split.col = paste0("fold", i),
              # alpha = 0.5,
              configs = cv_configs)$metric.val

outfile <- file.path(cv_configs$results.dir, 
                     paste0(weights_algo, "_screen",prop_keep, "_flat_metric.val.txt"))
write(out, file = outfile, ncolumns = 1)

cat("the cv results for inner fold ",i,"is stored in \n ",outfile,"\n")


