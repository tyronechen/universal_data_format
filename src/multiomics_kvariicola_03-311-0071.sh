Rscript run_pipeline.R \
   --classes ../data/sepsis/Kvariicola_03-311-0071/targets.tsv \
   --corr_cutoff 0.95 \
   --cross_val Mfold \
   --cross_val_folds 6 \
   --cross_val_nrepeat 24 \
   --dropna_classes FALSE \
   --dropna_prop 0 \
   --data ../data/sepsis/Kvariicola_03-311-0071/Metabolomics.tsv \
      ../data/sepsis/Kvariicola_03-311-0071/Proteomics.tsv \
      ../data/sepsis/Kvariicola_03-311-0071/Transcriptomics.tsv \
   --data_names metabolome proteome transcriptome \
   --force_unique FALSE \
   --ncpus 6 \
   --diablocomp 4 \
   --linkage 0.1 \
   --diablo_keepx 5 10 15 30 \
   --icomp 0 \
   --pcomp 10 \
   --plsdacomp 4 \
   --splsdacomp 4 \
   --splsda_keepx 5 10 15 30 \
   --dist_plsda mahalanobis.dist \
   --dist_splsda mahalanobis.dist \
   --dist_diablo mahalanobis.dist \
   --contrib max \
   --outfile_dir ../results/Kvariicola_03-311-0071/ \
   --rdata RData.RData \
   --plot Rplots.pdf \
   --args Rscript.sh
