#!/bin/bash
# download files required for setup
INFILE_DIR="../../data/rnaseq"

# download data from site http://bioinf.wehi.edu.au/RNAseqCaseStudy/
mkdir -p ${INFILE_DIR}
wget "http://bioinf.wehi.edu.au/RNAseqCaseStudy/data.tar.gz"
mv data.tar.gz ${INFILE_DIR}
cd ${INFILE_DIR}
tar -xzvf "data.tar.gz"
cd "../../src/"

# download the reference genome version we are using and reformat
# download gtf file
wget "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz"
mv "GCF_000001405.25_GRCh37.p13_genomic.gtf.gz" "${INFILE_DIR}"
# for this example we only use chromosome 1
# our version is slightly modified from the original gtf version
gzip -cd "${INFILE_DIR}/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz" | \
  grep -P "^#|^NC_000001.10\t" | sed "s|NC_000001.10	|chr1	|" > \
  "${INFILE_DIR}/hg19_chr1.gtf"

# download fasta file
# hg19_chr1.fa corresponds to this file
# wget 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz'
# gzip -d GCF_000001405.25_GRCh37.p13_genomic.fna.gz

echo "To generate the input bam files, follow the steps in http://bioinf.wehi.edu.au/RNAseqCaseStudy"
