# To run this case study, you should have R version of 3.0.2 or later 
# This is a modified version of the code available at:
#   http://bioinf.wehi.edu.au/RNAseqCaseStudy/
# Download and extract the files
# Running the 1.download.sh script should perform this automatically
# Then copy this script into that directory and run it instead of code.R
# load libraries
library(Rsubread)
library(limma)
library(edgeR)

# read in target file
options(digits=2)
targets <- readTargets()
targets

# create a design matrix
celltype <- factor(targets$CellType)
design <- model.matrix(~celltype)

# build an index for reference sequence (Chr1 in hg19)
buildindex(basename="chr1",reference="hg19_chr1.fa")

# align reads
align(index="chr1",readfile1=targets$InputFile,input_format="gzFASTQ",output_format="BAM",output_file=targets$OutputFile,unique=TRUE,indels=5)

# count numbers of reads mapped to NCBI Refseq genes
# the file is obtained from:
#   https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz
# our version is slightly modified from the original gtf version
# gzip -cd "${INFILE_DIR}/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz" | \
#   grep -P "^#|^NC_000001.10\t" | sed "s|NC_000001.10	|chr1	|" > \
#   "${INFILE_DIR}/hg19_chr1.gtf"

fc <- featureCounts(files=targets$OutputFile,isGTFAnnotationFile=TRUE,annot.ext="hg19_chr1.gtf")
write.table(fc$counts, sep="\t", quote=F, "original.tsv")
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

# generate RPKM values if you need them
x_rpkm <- rpkm(x,x$genes$Length)

# filter out low-count genes
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]

# perform voom normalization
y <- voom(x,design,plot=TRUE)

# cluster libraries
plotMDS(y,xlim=c(-2.5,2.5))

# fit linear model and assess differential expression
fit <- eBayes(lmFit(y,design))
top <- topTable(fit,coef=2,number=Inf,sort.by="P")
write.table(top, sep="\t", quote=F, file="original.tsv.top.tsv")
