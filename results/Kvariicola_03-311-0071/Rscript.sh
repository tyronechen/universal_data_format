Rscript run_pipeline.R \
   Rscript.sh  \
   --tune_off FALSE \
   --mini_run FALSE \
   --opts NA \
   --json NA \
   --classes ../data/sepsis/Kvariicola_03-311-0071/targets.tsv \
   --classes_secondary NA \
   --dropna_classes FALSE \
   --dropna_prop 0 \
   --data ../data/sepsis/Kvariicola_03-311-0071/Metabolomics.tsv ../data/sepsis/Kvariicola_03-311-0071/Proteomics.tsv ../data/sepsis/Kvariicola_03-311-0071/Transcriptomics.tsv \
   --data_names metabolome proteome transcriptome \
   --force_unique FALSE \
   --mappings NA \
   --ncpus 6 \
   --diablocomp 4 \
   --linkage 0.1 \
   --diablo_keepx 5 10 15 30 \
   --icomp 0 \
   --zero_as_na TRUE \
   --replace_missing TRUE \
   --pcomp 10 \
   --plsdacomp 4 \
   --splsdacomp 4 \
   --splsda_keepx 5 10 15 30 \
   --dist_plsda mahalanobis.dist \
   --dist_splsda mahalanobis.dist \
   --dist_diablo mahalanobis.dist \
   --cross_val Mfold \
   --cross_val_nrepeat 24 \
   --cross_val_folds 6 \
   --contrib max \
   --corr_cutoff 0.95 \
   --outfile_dir ../results/Kvariicola_03-311-0071/ \
   --rdata RData.RData \
   --plot Rplots.pdf
