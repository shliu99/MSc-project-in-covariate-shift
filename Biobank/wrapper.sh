# Bash wrapper

PHENO=breast
PROPBLACK=0.2
PROPSCREEN=0.05

for SDR in spca # kspca lfda klfda snmf
do
    for DRE in kmm #kliep ulsif
    do
        for FOLD in 1 2 3 4 5
        do
            qsub -N imp_${SDR}_${DRE}_${FOLD}_${PHENO}_${PROPBLACK}_${PROPSCREEN} \
              -v "prop_black=${PROPBLACK}" -v "sdr=${SDR}" -v "dre=${DRE}" \
              -v "pheno=${PHENO}" -v "fold=${FOLD}" -v "prop_screen=${PROPSCREEN}" \
              scripts/x04_importance_estimation_modified.R

            qsub -hold_jid imp_${SDR}_${DRE}_${FOLD}_${PHENO}_${PROPBLACK}_${PROPSCREEN} \
              -N init_${SDR}_${DRE}_${FOLD}_${PHENO}_${PROPBLACK}_${PROPSCREEN} \
              -v "prop_black=${PROPBLACK}" -v "sdr=${SDR}" -v "dre=${DRE}" \
              -v "pheno=${PHENO}" -v "fold=${FOLD}" -v "prop_screen=${PROPSCREEN}" \
              scripts/x05a_lasso_init_modified.R

            qsub -hold_jid init_${SDR}_${DRE}_${FOLD}_${PHENO}_${PROPBLACK}_${PROPSCREEN} \
              -N cv_${SDR}_${DRE}_${FOLD}_${PHENO}_${PROPBLACK}_${PROPSCREEN} \
              -v "prop_black=${PROPBLACK}" -v "sdr=${SDR}" -v "dre=${DRE}" \
              -v "pheno=${PHENO}" -v "fold=${FOLD}" -v "prop_screen=${PROPSCREEN}" \
              -t 5 -tc 1 scripts/x05b_lasso_cv_modified.R
            
            qsub -hold_jid cv_${SDR}_${DRE}_${FOLD}_${PHENO}_${PROPBLACK}_${PROPSCREEN} \
              -N fit_${SDR}_${DRE}_${FOLD}_${PHENO}_${PROPBLACK}_${PROPSCREEN} \
              -v "prop_black=${PROPBLACK}" -v "sdr=${SDR}" -v "dre=${DRE}" \
              -v "pheno=${PHENO}" -v "fold=${FOLD}" -v "prop_screen=${PROPSCREEN}" \
              scripts/x05c_lasso_final_modified.R
        done
    done
done