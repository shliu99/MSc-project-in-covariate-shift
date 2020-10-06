#!/apps/eb/skylake/software/R/3.6.2-foss-2019b/bin/Rscript
#$ -cwd
#$ -q short.qe
#$ -P mcvean.prjc
#$ -pe shmem 6
#$ -o logs/x05a_init_lasso.o.log
#$ -e logs/x05a_init_lasso.e.log
# This script initialises the lasso cross-validation for snpnet.

### Set up parameters ---

plink2_path  <- "/apps/eb/skylake/software/PLINK/2.00a2.3_x86_64/plink2"
zstdcat_path <- "/apps/eb/skylake/software/zstd/1.4.4-GCCcore-8.3.0/bin/zstdcat"

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
sdr <- Sys.getenv("sdr")
dre <- Sys.getenv("dre")
weights_algo = paste(sdr, dre, sep = "_")
gamma = 0.5
# 
# ncores <- as.integer(Sys.getenv("NSLOTS"))
# pheno  <- "breast"
# prop   <- 0.2
# fold   <- 1
# nfolds <- 5
# prop_keep = 0.05
# 
# sdr = "spca" # spca kspca lfda klfda snmf
# dre = "kmm" # kmm kliep ulsif
# weights_algo = paste(sdr, dre, sep = "_")

kfile  <- file.path("data", pheno,
                    paste0("fold", fold),
                    paste0("prop", prop, ".txt"))
pfile  <- "data/all_vars.tsv"
gfile  <- "data/ukb_cal_pgen/ukb_cal_all"
covars <- c("Sex", "Age", paste0("PC", 1:10), paste0("PC", 1:10, "_Sex"))

keep_ids <- read.table(kfile)[, 1]
ntrain <- length(keep_ids)

ids <- readIDsFromPsam(paste0(gfile, ".psam"))
id_df <- data.frame(ID = ids, stringsAsFactors = F) %>%
  mutate(sort_order = 1:n())

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
configs <- list(results.dir = file.path("temp", pheno,
                                        paste0("fold", fold),
                                        paste0("prop", prop)),
                plink2.path = plink2_path,
                zstdcat.path = zstdcat_path,
                keep = kfile)


for (f in seq_len(nfolds)) {
  dir.create(configs$results.dir,
             showWarnings = FALSE, recursive = TRUE)
}


configs$cv.full.prefix = file.path(configs$results.dir,
                                   paste0(weights_algo,
                                          "_screen",prop_keep,"_flat"))
### Prepare phenotype data.table ---

phe <- readPheMaster(phenotype.file = pfile, psam.ids = ids,
                     family = NULL, covariates = covars,
                     phenotype = pheno, status = NULL,
                     split.col = NULL, configs = configs)

folds <- cut(seq(1, nrow(phe)), breaks = nfolds, labels = FALSE)
foldids <- sample(folds, length(folds))

cv.phenotype.file <- paste0(configs[["cv.full.prefix"]], ".tsv")
phe %>%
  dplyr::mutate(fold = foldids, tmp = "val") %>%
  tidyr::pivot_wider(names_prefix = "fold", names_from = fold,
                     values_from = tmp, values_fill = list(tmp = "train")) %>%
  dplyr::mutate(weights = weights) %>%
  data.table::fwrite(cv.phenotype.file, sep = "\t")

### Pre-compute lambda path ---

gc()

if (pheno %in% c("breast", "cancer")) {
  covars <- c("Age", paste0("PC", 1:10))
} else {
  covars <- c("Sex", "Age", paste0("PC", 1:10), paste0("PC", 1:10, "_Sex"))
}

configs$covariates <- covars
configs$nCores <- ncores
full.lams <- snpnet(genotype.pfile = gfile,
                    phenotype.file = pfile,
                    phenotype = pheno,
                    covariates = covars,
                    weights = weights,
                    configs = configs,
                    lambda_only = TRUE)

saveRDS(configs, file.path(configs$results.dir, 
                           paste0(weights_algo, "_screen",prop_keep,"_flat_configs.RDS")))

full.lams.file <- file.path(configs$results.dir, paste0(weights_algo,
                                                        "_screen",prop_keep,
                                                        "_flat_full.lams.txt"))
write(full.lams, file = full.lams.file, ncolumns = 1)
