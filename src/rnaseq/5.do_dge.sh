#!/bin/bash
for i in ../../results/rnaseq/*joined.tsv;
do
  outfile_dir=$(dirname ${i})
  Rscript ../differential_gene_expression.r \
    ${i} \
    ../../data/rnaseq/Targets.txt \
    -o ${outfile_dir}/ \
    -n ${i/joined/dge}
done
