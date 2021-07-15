#!/usr/bin/Rscript
library(ggplot2)
library(Rsubread)
library(limma)
library(edgeR)

# read in target file
options(digits=2)
targets <- readTargets("../../data/rnaseq/Targets.txt")

# create a design matrix
celltype <- factor(targets$CellType)
design <- model.matrix(~celltype)

# build an index for reference sequence (Chr1 in hg19)
# buildindex(basename="chr1",reference="../../data/rnaseq/hg19_chr1.fa")
#
# # align reads
# align(
#   ndex="chr1",readfile1=targets$InputFile,input_format="gzFASTQ",
#   output_format="BAM",output_file=targets$OutputFile,unique=TRUE,indels=5
# )

# count numbers of reads mapped to NCBI Refseq genes
# the file is obtained from:
#   https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz
# our version is slightly modified from the original gtf version
# gzip -cd "${INFILE_DIR}/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz" | \
#   grep -P "^#|^NC_000001.10\t" | sed "s|NC_000001.10	|chr1	|" > \
#   "${INFILE_DIR}/hg19_chr1.gtf"
fc <- featureCounts(
  files=paste(
    "../../data/rnaseq/", c("A_1.bam", "A_2.bam", "B_1.bam", "B_2.bam"), sep=""
  ),
  isGTFAnnotationFile=TRUE,
  annot.ext="../../data/rnaseq/hg19_chr1.gtf.gz"
)
# write.table(fc$counts, sep="\t", quote=F, "original.tsv")
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

# filter out low-count genes
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]

# generate CPM values for plotting
x_cpm <- cpm(x)

# make pca plots
y <- prcomp(t(x_cpm), scale=TRUE)
df_out <- as.data.frame(y$x)
df_out$group <- sapply(
  strsplit(as.character(row.names(df_out)), ".", fixed=TRUE), "[[", 8
)
df_out$group <- sapply( strsplit(df_out$group, "_", fixed=TRUE), "[[", 1 )

p <- ggplot(df_out, aes(x=PC1, y=PC2, color=group)) + geom_point(size=8) + theme_bw()
ggsave("../../results/ismb_plots_2021/original.pca.pdf", p)

# clear environment to avoid variable clash
rm(list=ls())

# read in target file
targets <- readTargets("../../data/rnaseq/Targets.txt")

# create a design matrix
celltype <- factor(targets$CellType)
design <- model.matrix(~celltype)

# make pca plot from reformatted data
# data was reformatted by passing it through our pipeline
x <- read.table("../../results/rnaseq/24000.joined.tsv")
x <- DGEList(x)

# filter out low-count genes
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]

# generate CPM values for plotting
x_cpm <- cpm(x)

# make pca plots
y <- prcomp(t(x_cpm), scale=TRUE)
df_out <- as.data.frame(y$x)
df_out$group <- sapply( strsplit(row.names(df_out), "_", fixed=TRUE), "[[", 1 )

p <- ggplot(df_out, aes(x=PC1, y=PC2, color=group)) + geom_point(size=8) + theme_bw()
ggsave("../../results/ismb_plots_2021/reformat.pca.pdf", p)
