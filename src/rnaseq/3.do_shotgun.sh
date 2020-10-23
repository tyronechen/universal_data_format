#!/bin/bash
genome_base="hg19_chr1"
outfile_dir="../../results/rnaseq/"

for i in ../../data/rnaseq/*bam
do
  j=$(basename $i)
  for k in 100 1000 10000 24000 100000
  do
    ../shotgun_genome.sh \
      ../../data/rnaseq/${genome_base}.fa \
      ${k} \
      ${outfile_dir}/${j}.${k}
    ../get_coverage.sh \
      ${outfile_dir}/${j}.${k}/${genome_base}.bed \
      ../../data/rnaseq/${genome_base}.fa.txt \
      ${i} \
      ${outfile_dir}/${j}.${k}
    #Rscript ../differential_gene_expression.r \
    #  ${outfile_dir}/${j}.${k}/${j}.cov.counts.tsv \
    #  ../../data/rnaseq/Targets.txt \
    #  -o ${outfile_dir}/${j}.${k}/ \
    #  -n ${j}.${k}
  done
done
