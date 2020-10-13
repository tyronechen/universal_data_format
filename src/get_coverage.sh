#!/bin/bash
# take reference genome, take data bed files, return counts per genomic region
#   $1 = genome_path
#   $2 = genome_size
#   $3 = sample_name
#   $4 = outfile_dir
if [[ $1 == "" ]] || [[ $2 == "" ]]; then
  echo "Usage: get_coverage <GENOME_PATH> <GENOME_SIZE> <SAMPLE_PATH> [OUTFILE_DIR]" | fmt
  echo "Genomes should be available as bed, genome size as txt, samples as bed." | fmt
  exit 1
else
  GENOME_SIZE=${2}
  BASE_GENOME_SIZE=$(basename ${2})
  FULL_GENOME_PATH=${1}
  FULL_SAMPLE_PATH=${3}
  BASE_GENOME_PATH=$(basename ${1})
  BASE_SAMPLE_PATH=$(basename ${3})
  GENOME_PATH=$(basename ${1} | cut -d . -f1)
  SAMPLE_PATH=$(basename ${3} | cut -d . -f1)
  echo "# Genome name:" ${GENOME_PATH}
  echo "# Genome size:" ${GENOME_SIZE}
  echo "# Sample name:" ${SAMPLE_PATH}
fi

if [[ $4 == "" ]] || [[ $4 == "." ]] || [[ $4 == "./" ]]; then
  OUTFILE_DIR="./out/"
else
  OUTFILE_DIR=${4}
fi
echo "# Outfile dir:" ${OUTFILE_DIR}
mkdir -p ${OUTFILE_DIR}

check_software() {
  # install necessary software utilities from:
  # http://hgdownload.soe.ucsc.edu/admin/exe/ and bedtools
  which bamToBed > /dev/null
  if [[ $? == "" ]]; then
    echo "# Download bamToBed from http://hgdownload.soe.ucsc.edu/admin/exe/"
    exit 1
  else
    echo "# Calling bamToBed from this path:"
    echo "#" $(which bamToBed)
  fi
  # install bedtools from source or conda install bedtools
  which bedtools > /dev/null
  if [[ $? == "" ]]; then
    echo "# Download bedtools from source or conda install bedtools"
    exit 1
  else
    echo "# Calling bedtools from this path:"
    echo "#" $(which bamToBed)
  fi
}

bam_to_bed() {
  # transform bamfile into bedfile
  #   bam_to_bed /path/to/bam /path/to/bed: -> bed unsorted
  echo "bamToBed -i ${1}" '>' "${2}"
  bamToBed -i ${1} > ${2}
}

sort_bed() {
  # sort bedfile data by chromosome index
  #   sort_bed /path/to/bed /path/to/genome.txt /path/to/out: -> bed
  # add this if want bytesort instead of numeric: sort -T ./ -k1,1 -k2,2n
  echo "bedtools sort -i ${1} -g ${2}" '>' "${3}"
  bedtools sort -i ${1} -g ${2} > ${3}
  if [[ $? == 0 ]]; then rm ${1}; fi
}

get_coverage() {
  # obtain genome coverage proportion (normalises within modality)
  #   get_coverage /path/to/data.bed /path/to/genome.bed /path/to/genome.txt /path/to/out: -> coverage
  # NOTE: if unsorted, you will probably get a memory leak!
  echo "bedtools coverage -sorted -g ${3} -a ${2} -b ${1}" '>' "${4}"
  bedtools coverage -sorted -g ${3} -a ${2} -b ${1} > ${4}
  if [[ $? == 0 ]]; then rm ${1}; fi
}

write_counts() {
  # write counts per block out to file
  #   write_counts /path/to/data.cov: -> counts
  echo 'printf' "\t${SAMPLE_PATH}\n" '>' "${1}.counts.tsv"
  echo "paste <(cut -f1-3 ${1} | tr '\t' '_') <(cut -f5 ${1})" '>>' "${1}.counts.tsv"
  printf "\t${SAMPLE_PATH}\n" > ${1}.counts.tsv
  paste <(cut -f1-3 ${1} | tr '\t' '_') <(cut -f5 ${1}) >> ${1}.counts.tsv
  if [[ $? == 0 ]]; then rm ${1}; fi
}

main() {
  check_software
  # data operations
  # bam_to_bed BAMFILE OUTFILE_PATH: -> UNSORTED BED
  bam_to_bed ${FULL_SAMPLE_PATH} \
    ${OUTFILE_DIR}/${BASE_SAMPLE_PATH}.bed.unsort

  # sort_bed UNSORTED.BED GENOME.TXT: OUTFILE_PATH -> SORTED.BED
  sort_bed ${OUTFILE_DIR}/${BASE_SAMPLE_PATH}.bed.unsort \
    ${GENOME_SIZE} \
    ${OUTFILE_DIR}/${BASE_SAMPLE_PATH}.bed

  # get_coverage SORTED.BED GENOME.BED OUTFILE_PATH: -> COVERAGE.COV
  get_coverage ${OUTFILE_DIR}/${BASE_SAMPLE_PATH}.bed \
    ${FULL_GENOME_PATH/%.bed}.fa.bed \
    ${GENOME_SIZE} \
    ${OUTFILE_DIR}/${BASE_SAMPLE_PATH}.cov
  
  # get_bedgraph ${SAMPLE_PATH}.bed ${GENOME_PATH}
  write_counts ${OUTFILE_DIR}/${BASE_SAMPLE_PATH}.cov
}

main
