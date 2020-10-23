#!/bin/bash
INFILE_DIR="../../results/rnaseq/"
for i in 100 1000 10000 24000 100000
do
  files=$(find ${INFILE_DIR}/*$i/ -name "*cov.counts.tsv" -type f)
  python ../join_counts.py $files -o ${INFILE_DIR}/${i}.joined.tsv
done
