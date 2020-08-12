#!/bin/bash
# take reference genome, take data bed files, return normalised genome coverage
#   $1 = genome_name
#   $2 = sample_name
if [[ $1 == "" ]] || [[ $2 == "" ]]; then
  echo "Usage: get_coverage <genome_name> <sample_name>" | fmt
  echo "Genomes should be available as txt and bed, samples as bed." | fmt
  exit 1
else
  GENOME_NAME=$(basename ${1} | cut -d . -f1)
  SAMPLE_NAME=$(basename ${2} | cut -d . -f1)
  echo "Genome name:" ${GENOME_NAME}
  echo "Sample name:" ${SAMPLE_NAME}
fi

check_software() {
  # install necessary software utilities from:
  # http://hgdownload.soe.ucsc.edu/admin/exe/ and bedtools
  which bamToBed
  if [[ $? == "" ]]; then
    echo "Download bamToBed from http://hgdownload.soe.ucsc.edu/admin/exe/"
    exit 1
  else
    echo "Calling bamToBed from this path:"
  fi
  # install bedtools from source or conda install bedtools
  which bedtools
  if [[ $? == "" ]]; then
    echo "Download bedtools from source or conda install bedtools"
    exit 1
  else
    echo "Calling bedtools from this path:"
  fi
}

get_genome() {
  # download matching genome assembly from link
  # NOTE: remember to check bam and bedfile for custom chromosome identifiers
  #   get_genome /url/to/genome: -> genome.fa
  genome=$(echo ${1} | sed "s|.*/||")
  wget ${1}
  gzip -dv ${genome}
}

genome_to_txt() {
  # convert genome to bedtools compatible format
  #   genome_to_bed /path/to/genome.fa: -> index
  faidx -i bed ${1} | cut -f1,3 | grep -v "_" | grep -v "EBV" > ${1/.fasta/}.txt
}

shotgun_genome() {
  # shotgun_genome genome into N length blocks to get higher resolution coverage
  #   shotgun_genome /path/to/genome/index.txt 1000: -> bed
  echo bedtools makewindows -g ${1} -w ${2} -i winnum ${1/.txt/}.bed
  bedtools makewindows -g ${1} -w ${2} -i winnum > ${1/.txt/}.bed
}

bam_to_bed() {
  # transform bamfile into bedfile
  #   bam_to_bed /path/to/bam: -> bed unsorted
  echo bamToBed -i ${1} ${1/.bam/}.bed.unsort
  bamToBed -i ${1} > ${1/.bam/}.bed.unsort
}

sort_bed() {
  # sort bedfile data by chromosome index
  #   sort_bed /path/to/bed /path/to/genome.txt: -> bed
  # add this if want bytesort instead of numeric: sort -T ./ -k1,1 -k2,2n
  echo bedtools sort -i ${1} -g ${2} ${1/.unsort/}
  bedtools sort -i ${1} -g ${2} > ${1/.unsort/}
}

get_coverage() {
  # obtain genome coverage proportion (normalises within modality)
  #   get_coverage /path/to/data.bed /path/to/genome_name: -> coverage
  # NOTE: if unsorted, you will probably get a memory leak!
  echo bedtools coverage -sorted -g ${2}.txt -a ${2}.bed -b ${1} ${1/.bed/}.cov
  bedtools coverage -sorted -g ${2}.txt -a ${2}.bed -b ${1} > ${1/.bed/}.cov
}

filter_zero() {
  # filter regions empty across all samples?
  # get last col for each coverage set and remove?
  echo ""
}

main() {
  check_software
  # genome operations
  # get_genome URL: -> GENOME.FASTA
  # genome_to_txt GENOME.FASTA: -> INDEX.TXT
  #   genome_to_txt hg19.fasta
  # shotgun_genome INDEX.TXT WINDOW: -> GENOME.BED
  #   shotgun_genome hg19.txt 1000

  # data operations
  # bam_to_bed BAMFILE: -> UNSORTED BED
  # sort_bed UNSORTED.BED GENOME.TXT: -> SORTED.BED
  # sort_bed ${SAMPLE_NAME}.bed.unsort ${GENOME_NAME}.txt
  # get_coverage SORTED.BED GENOME.BED: -> COVERAGE.COV
  # get_coverage ${SAMPLE_NAME}.bed ${GENOME_NAME}
  get_bedgraph ${SAMPLE_NAME}.bed ${GENOME_NAME}
}

main
