#!/bin/bash
# take reference genome, take data bed files, return counts per genomic region
#   $1 = genome_name
#   $2 = sample_name
#   $3 = window_size
if [[ $1 == "" ]] || [[ $2 == "" ]]; then
  echo "Usage: get_coverage <GENOME_SIZE> <SAMPLE_NAME> [WINDOW_SIZE]" | fmt
  echo "Genomes should be available as txt and bed, samples as bed." | fmt
  exit 1
else
  FULL_GENOME_NAME=${1}
  FULL_SAMPLE_NAME=${2}
  BASE_GENOME_NAME=$(basename ${1})
  BASE_SAMPLE_NAME=$(basename ${2})
  GENOME_NAME=$(basename ${1} | cut -d . -f1)
  SAMPLE_NAME=$(basename ${2} | cut -d . -f1)
  echo "# Genome name:" ${GENOME_NAME}
  echo "# Sample name:" ${SAMPLE_NAME}
fi

if [[ $3 == "" ]]; then
  WINDOW_SIZE=1000
else
  WINDOW_SIZE=${3}
fi
echo "# Window size:" ${WINDOW_SIZE} "(default: 1000)"

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
  # install samtools from source or conda install samtools
  which samtools > /dev/null
  if [[ $? == "" ]]; then
    echo "# Download samtools from source or conda install samtools"
    exit 1
  else
    echo "# Calling samtools from this path:"
    echo "#" $(which samtools)
  fi
}

get_genome() {
  # download matching genome assembly from link
  # NOTE: remember to check bam and bedfile for custom chromosome identifiers
  #   get_genome /url/to/genome: -> genome.fa
  genome=$(echo ${1} | sed "s|.*/||")
  echo wget ${1}
  wget ${1}
  echo gzip -dv ${genome}
  gzip -dv ${genome}
}

genome_to_txt() {
  # convert genome to bedtools compatible format
  #   genome_to_bed /path/to/genome.fa: -> index
  echo "faidx -i bed ${1} | cut -f1,3 | grep -v "_" | grep -v EBV" '>' "${1}.txt"
  faidx -i bed ${1} | cut -f1,3 | grep -v "_" | grep -v EBV > ${1}.txt
}

shotgun_genome() {
  # shotgun_genome genome into N length blocks to get higher resolution coverage
  #   shotgun_genome /path/to/genome/index.txt 1000 outfile: -> bed
  echo "bedtools makewindows -g ${1} -w ${2} -i winnum" '>' "${3}"
  bedtools makewindows -g ${1} -w ${2} -i winnum > ${3}
}

bam_to_bed() {
  # transform bamfile into bedfile
  #   bam_to_bed /path/to/bam: -> bed unsorted
  outfile_path=$(basename ${1} | cut -d . -f1)
  echo "bamToBed -i ${1}" '>' "${outfile_path}.bed.unsort"
  bamToBed -i ${1} > ${outfile_path}.bed.unsort
}

sort_bed() {
  # sort bedfile data by chromosome index
  #   sort_bed /path/to/bed /path/to/genome.txt: -> bed
  # add this if want bytesort instead of numeric: sort -T ./ -k1,1 -k2,2n
  echo "bedtools sort -i ${1} -g ${2}" '>' "${1/.unsort/}"
  bedtools sort -i ${1} -g ${2} > ${1/.unsort/}
}

get_coverage() {
  # obtain genome coverage proportion (normalises within modality)
  #   get_coverage /path/to/data.bed /path/to/genome.bed /path/to/genome.txt: -> coverage
  # NOTE: if unsorted, you will probably get a memory leak!
  echo "bedtools coverage -sorted -g ${2} -a ${2} -b ${1}" '>' "${1/.bed/}.cov"
  bedtools coverage -sorted -g ${3} -a ${2} -b ${1} > ${1/.bed/}.cov
}

write_counts() {
  # write counts per block out to file
  #   write_counts /path/to/data.cov: -> counts
  echo 'printf' "Features\t${1/.cov/}\n" '>' "${1}.counts.tsv"
  echo "paste <(cut -f1-3 ${1} | tr '\t' '_') <(cut -f5 ${1})" '>>' "${1}.counts.tsv"
  printf "Features\t${1/.cov/}\n" > ${1}.counts.tsv
  paste <(cut -f1-3 ${1} | tr '\t' '_') <(cut -f5 ${1}) >> ${1}.counts.tsv
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
  genome_to_txt ${FULL_GENOME_NAME} # hg19_chr1.fa
  # shotgun_genome INDEX.TXT WINDOW: -> GENOME.BED
  shotgun_genome ${FULL_GENOME_NAME}.txt ${WINDOW_SIZE} ${BASE_GENOME_NAME}.bed

  # data operations
  # bam_to_bed BAMFILE: -> UNSORTED BED
  bam_to_bed ${FULL_SAMPLE_NAME}
  # sort_bed UNSORTED.BED GENOME.TXT: -> SORTED.BED
  sort_bed ${SAMPLE_NAME}.bed.unsort ${FULL_GENOME_NAME}.txt
  # get_coverage SORTED.BED GENOME.BED: -> COVERAGE.COV
  get_coverage ${SAMPLE_NAME}.bed ${BASE_GENOME_NAME}.bed ${FULL_GENOME_NAME}.txt
  # get_bedgraph ${SAMPLE_NAME}.bed ${GENOME_NAME}
  # write_counts ${SAMPLE_NAME}.cov
  write_counts ${SAMPLE_NAME}.cov
  
  # cleanup intermediate files
  mv ${SAMPLE_NAME}.cov.counts.tsv ${OUTFILE_DIR}/${SAMPLE_NAME}.counts.tsv
  if [[ $? == 0 ]]; then
    rm ${SAMPLE_NAME}.bed.unsort \
      ${SAMPLE_NAME}.bed \
      ${SAMPLE_NAME}.cov \
      ${BASE_GENOME_NAME}.bed
  fi
}

main
