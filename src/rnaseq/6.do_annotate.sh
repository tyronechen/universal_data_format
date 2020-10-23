#!/bin/bash
INFILE_DIR="../../results/rnaseq/"
for i in 100 1000 10000 24000 100000
do
  python ../annotate_regions.py "${INFILE_DIR}/${i}.dge.tsv.top.tsv" \
    "../../data/rnaseq/hg19_chr1.gtf.gz" \
    -o "${INFILE_DIR}/${i}.annotated.tsv"
done
