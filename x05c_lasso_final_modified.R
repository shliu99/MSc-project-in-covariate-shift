#!/apps/eb/skylake/software/R/3.6.2-foss-2019b/bin/Rscript
#$ -cwd
#$ -q short.qe
#$ -P mcvean.prjc
#$ -m ea
#$ -pe shmem 4
#$ -o logs/x05c_lasso_final.o.log
#$ -e logs/x05c_lasso_final.e.log

### Set up parameters ---

# devtools::install_github("brieuclehmann/snpnet")
library(snpnet)
library(readr)
library(dplyr)
set.seed(1)

ncores <- as.integer(Sys.getenv("NSLOTS"))
pheno <- Sys.getenv("pheno")
prop <- as.double(Sys.getenv("prop_black"))
prop_keep <- as.double(Sys.getenv("prop_screen"))
fold <- as.integer(Sys.getenv("fold"))
nfolds <- 5
i <- as.integer(Sys.getenv("SGE_TASK_ID"))
sdr <- Sys.getenv("sdr")
dre <- Sys.getenv("dre")
weights_algo = paste(sdr, dre, sep = "_")
gamma = 0.5

# ################################################################################
# ncores <- as.integer(Sys.getenv("NSLOTS"))
# pheno  <- "breast"
# prop   <- 0.2
# fold   <- 1
# nfolds <- 5
# prop_keep = 0.05
# sdr = "default" # spca kspca lfda klfda snmf
# dre = "default" # kmm kliep ulsif
# weights_algo = paste(sdr, dre, sep = "_")
# ################################################################################

kfile  <- file.path("data", pheno,
                    paste0("fold", fold),
                    paste0("prop", prop, ".txt"))

pfile  <- "data/all_vars.tsv"
gfile  <- "data/ukb_cal_pgen/ukb_cal_all"

results.dir <- file.path("temp", pheno,
                         paste0("fold", fold),
                         paste0("prop", prop))
configs   <- readRDS(file.path(results.dir, 
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
weights <- as.numeric(readLines(file.path(weights_path,paste0(weights_algo,
                                                              "_screen",prop_keep,
                                                              ".txt"))))
weights <- weights^gamma
weights = weights/mean(weights)
weights[weights > 5] <- 5
weights[weights <= 1e-3] <- 0
################################################################################
### Gather cross-validation results ---

inner_folds <- seq_len(nfolds)
cvout <- inner_folds %>%
  sapply(function(x) file.path(results.dir, paste0(paste0("foldi", x)),
                               paste0(weights_algo, "_screen",prop_keep, "_flat_metric.val.txt"))) %>%
  sapply(scan)

cvm <- apply(cvout, 1, mean)
cvsd <- apply(cvout, 1, sd)
lambda.min <- full.lams[which.max(cvm)]
lambda_na <- apply(cvout, 1, function(x) !all(is.na(x)))

### Fit lasso on full dataset ---

fit.lams <- full.lams[lambda_na]
configs$nCores <- ncores

snpnet.object <- snpnet(gfile,
                        pfile,
                        pheno,
                        covariates = covars,
                        weights = weights,
                        full.lams = fit.lams,
                        configs = configs,
                        mem = 8e3 * ncores)

if (is.na(cvm[which.max(cvm) + 1]))
  warning("Cross-validation may have stopped early. Consider increasing stopping.lag.")

### Save output ---

mod <- list(snpnet.object = snpnet.object,
            cvm = cvm,
            cvsd = cvsd,
            cvout = cvout,
            lambda.min = lambda.min,
            full.lams = full.lams,
            fit.lams = fit.lams)

outdir <- file.path("output", pheno, paste0("fold", fold), paste0("prop", prop))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, paste0(weights_algo,
                                    "_screen",prop_keep,
                                    "_flat.RDS"))
saveRDS(mod, outfile)


### Predictions ---

# Prepare test set
train_id <- read.table(kfile)[, 1]
pheno_df <- read_tsv(pfile) %>%
  mutate(ID = paste(FID, IID, sep = "_"), train = (eid %in% train_id)) %>%
  rename(y = all_of(pheno))


kfile_test <- file.path("data", pheno,
                       paste0("fold", fold),
                       paste0("prop", prop, "_test.txt"))


gfile_test <- file.path("data", pheno,
                        paste0("fold", fold),
                        paste0("prop", prop, "_ukb_cal_test"))

cmd_plink2 <- paste("/apps/eb/skylake/software/PLINK/2.00a2.3_x86_64/plink2",
                    "--pfile", gfile,
                    "--threads", ncores,
                    "--keep", kfile_test,
                    "--make-pgen vzs",
                    "--out", gfile_test)
system(cmd_plink2, intern = F)


# Get predictions for best lambda
lambda_ind <- which(mod$full.lams == mod$lambda.min)

configs <- list(plink2.path = "/apps/eb/skylake/software/PLINK/2.00a2.3_x86_64/plink2",
                zstdcat.path = "/apps/eb/skylake/software/zstd/1.4.4-GCCcore-8.3.0/bin/zstdcat")

pred <- predict_snpnet(snpnet.object,
                       new_genotype_file = gfile_test,
                       new_phenotype_file = pfile,
                       phenotype = pheno,
                       covariate_names = covars,
                       idx = lambda_ind,
                       configs = configs)

pred_df <- tibble(ID = rownames(pred$prediction$train),
                  pred = pred$prediction$train[, 1]) %>%
  inner_join(pheno_df, by = "ID") %>%
  mutate(pred_bin = 1 / (1 + exp(-pred)), type = "lasso")

pred_file = file.path(outdir,paste0(weights_algo, "_screen",prop_keep,"_flat_res.tsv"))

write_tsv(pred_df, pred_file)

file.remove(list.files(dirname(gfile_test), basename(gfile_test),
                       full.names = TRUE))
