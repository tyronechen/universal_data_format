#!/bin/bash
# take reference genome, return genome divided into blocks
#   $1 = genome_path
#   $2 = window_size
#   $3 = outfile_dir
if [[ $1 == "" ]] ; then
  echo "Usage: shotgun_genome.sh <GENOME_PATH> [WINDOW_SIZE] [OUTFILE_DIR]" | fmt
  echo "Genome should be available as fasta file." | fmt
  exit 1
else
  FULL_GENOME_PATH=${1}
  BASE_GENOME_PATH=$(basename ${1})
  GENOME_PATH=$(basename ${1} | cut -d . -f1)
  GENOME_DIR=$(dirname ${FULL_GENOME_PATH})
  echo "# Genome name:" ${GENOME_PATH}
fi

if [[ $2 == "" ]]; then
  WINDOW_SIZE=1000
else
  WINDOW_SIZE=${2}
fi
echo "# Window size:" ${WINDOW_SIZE} "(default: 1000)"

if [[ $3 == "" ]] || [[ $3 == "." ]] || [[ $3 == "./" ]]; then
  OUTFILE_DIR="./out/"
else
  OUTFILE_DIR=${3}
fi
echo "# Outfile dir:" ${OUTFILE_DIR}
mkdir -p ${OUTFILE_DIR}

check_software() {
  # install necessary software utilities from:
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
  echo "faidx -i bed ${1} | cut -f1,3" '>' "${1}.txt" #${GENOME_DIR}/${BASE_GENOME_PATH}.txt"
  faidx -i bed ${1} | cut -f1,3 > ${1}.txt #${GENOME_DIR}/${BASE_GENOME_PATH}.txt
}

shotgun_genome() {
  # shotgun_genome genome into N length blocks to get higher resolution coverage
  #   shotgun_genome /path/to/genome/index.txt 1000: -> bed
  echo "bedtools makewindows -g ${1} -w ${2} -i winnum" '>' "${OUTFILE_DIR}/${BASE_GENOME_PATH}.bed"
  bedtools makewindows -g ${1} -w ${2} -i winnum > ${OUTFILE_DIR}/${BASE_GENOME_PATH}.bed
}

main() {
  check_software
  # genome operations
  # get_genome URL: -> GENOME.FASTA
  # genome_to_txt GENOME.FASTA: -> INDEX.TXT
  genome_to_txt ${FULL_GENOME_PATH} # hg19_chr1.fa
  # shotgun_genome INDEX.TXT WINDOW: -> GENOME.BED
  shotgun_genome ${FULL_GENOME_PATH}.txt ${WINDOW_SIZE} #${OUTFILE_DIR}/${BASE_GENOME_PATH}.txt ${WINDOW_SIZE}
}

main
