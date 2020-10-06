#!/apps/eb/skylake/software/R/3.6.2-foss-2019b/bin/Rscript
#$ -cwd
#$ -P mcvean.prjc
#$ -N prep_pheno
#$ -q short.qe
#$ -o logs/prep_pheno.log
#$ -e logs/prep_pheno.log


pheno <- "breast" # HC HT T2D asthma prostate breast
prop_min <- 0.2
nfolds <- 5

library(Matching)
library(readr)
library(dplyr)
set.seed(1) # for reproducibility

min_coding <- 4 # Black or Black British
maj_coding <- 1 # White

all_of <- function (x) {
  if (is.function(x)) {
    vctrs::vec_as_location(x, 0L)
    abort("Internal error: `all_of()` should have failed sooner")
  }
  x
}

pheno_df <- read_tsv("data/all_vars.tsv") %>%
  rename(y = all_of(pheno)) %>%
  filter(!is.na(y) & parent_coding %in% c(min_coding, maj_coding)) %>%
  mutate(Tr = (parent_coding == min_coding))

for (f in seq_len(nfolds)) {
  dir.create(file.path("data", pheno, paste0("fold", f)),
             showWarnings = FALSE, recursive = TRUE)
}

# Train-test split
train_ratio <- round((1 - prop_min) / prop_min)

min_df <- pheno_df %>%
  filter(parent_coding == min_coding) %>%
  mutate(fold = sample(cut(seq_len(n()),
                           breaks = nfolds, labels = FALSE))) %>%
  dplyr::select(eid, Tr, Sex, Age, fold)

maj_df <- pheno_df %>%
  filter(parent_coding == maj_coding) %>%
  mutate(fold = NA_integer_) %>%
  dplyr::select(eid, Tr, Sex, Age, fold)

for (f in seq_len(nfolds)) {
  # train files
  match_df <- min_df %>%
    filter(fold != f) %>%
    bind_rows(maj_df)

  match_out <- match_df %>%
    with(Match(Tr = Tr, X = cbind(Sex, Age), exact = c(TRUE, FALSE),
               M = train_ratio, Z = eid, ties = FALSE, replace = FALSE))

  idx <- unique(c(match_out$index.treated, match_out$index.control))
  train_id <- match_df$eid[idx]

  file_train <- file.path("data", pheno,
                       paste0("fold", f),
                       paste0("prop", prop_min, ".txt"))
  write(rep(train_id, each = 2), file_train, ncol = 2)

  # test files
  test_id <- min_df %>%
    filter(fold == f) %>%
    pull(eid)

  file_test <- file.path("data", pheno,
                          paste0("fold", f),
                          paste0("prop", prop_min, "_test.txt"))

  write(rep(test_id, each = 2), file_test, ncol = 2)
}

