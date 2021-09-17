Rscript run_pipeline.R \
   --classes ../data/sepsis/Ecoli_B36/targets.tsv \
   --corr_cutoff 0.95 \
   --cross_val Mfold \
   --cross_val_folds 2 \
   --cross_val_nrepeat 10 \
   --dropna_classes FALSE \
   --dropna_prop 0 \
   --data ../data/sepsis/Ecoli_B36/Metabolomics.tsv \
      ../data/sepsis/Ecoli_B36/Proteomics.tsv \
      ../data/sepsis/Ecoli_B36/Transcriptomics.tsv \
   --data_names metabolome proteome transcriptome \
   --force_unique FALSE \
   --ncpus 1 \
   --diablocomp 4 \
   --linkage 0.1 \
   --diablo_keepx 5 10 15 30 \
   --icomp 0 \
   --pcomp 10 \
   --plsdacomp 2 \
   --splsdacomp 3 \
   --splsda_keepx 5 10 15 \
   --dist_plsda mahalanobis.dist \
   --dist_splsda mahalanobis.dist \
   --dist_diablo mahalanobis.dist \
   --contrib max \
   --outfile_dir ../results/Ecoli_B36/ \
   --rdata RData.RData \
   --plot Rplots.pdf \
   --args Rscript.sh
