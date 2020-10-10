#!/usr/bin/Rscript
library(argparser)
library(edgeR)
library(Rsubread)

parse_argv <- function() {
    library(argparser, quietly=TRUE)
    p <- arg_parser("Run a differential expression analysis.")

    # Add command line arguments
    p <- add_argument(p, "infile", type="character", help="sample data")
    p <- add_argument(p, "targets", type="character", help="sample information")
    p <- add_argument(p, "--outfile_dir", type="character", default="./",
                      help="sample information")
    p <- add_argument(p, "--plot", type="character", help="output plots file")
    # Parse the command line arguments
    argv <- parse_args(p)

    # Do work based on the passed arguments
    return(argv)
}

main <- function() {
  argv <- parse_argv()
  pdf(paste(argv$outfile_dir, "/", argv$plot, sep=""))
  targets <- readTargets(argv$targets)

  # create a design matrix
  celltype <- factor(targets$CellType)
  design <- model.matrix(~celltype)

  # count numbers of reads mapped to NCBI Refseq genes
  data <- read.table(argv$infile, sep="\t")
  dge <- DGEList(counts=data)

  # filter out low-count genes
  isexpr <- rowSums(cpm(dge) > 10) >= 2
  dge <- dge[isexpr,]

  # perform voom normalization
  voomed <- voom(dge, design, plot=TRUE)

  # cluster libraries
  plotMDS(voomed, xlim=c(-2.5,2.5))

  # fit linear model and assess differential expression
  fit <- eBayes(lmFit(voomed, design))
  topTable(fit, coef=2)

  dev.off()
}

main()
