#!/usr/bin/Rscript
# This is a modified version of the code found here:
#   http://bioinf.wehi.edu.au/RNAseqCaseStudy/

# To run this case study, you should have R version of 3.0.2 or later
library(argparser)
library(edgeR)
library(Rsubread)

parse_argv <- function() {
    p <- argparser::arg_parser("Run a differential expression analysis.")

    # Add command line arguments
    p <- argparser::add_argument(p, "infile", type="character",
                                 help="sample data")
    p <- argparser::add_argument(p, "targets", type="character",
                                 help="sample information")
    p <- argparser::add_argument(p, "--outfile_dir", type="character",
                                 default="./", help="output file directory")
    p <- argparser::add_argument(p, "--names", type="character",
                                 help="output names")
    p <- argparser::add_argument(p, "--voom_off", flag=TRUE,
                                 help="toggle voom off")
    p <- argparser::add_argument(p, "--ebayes_off", flag=TRUE,
                                 help="toggle empirical bayes off")
    # Parse the command line arguments
    argv <- argparser::parse_args(p)

    # Do work based on the passed arguments
    return(argv)
}

write_args <- function(args, argpath) {
  if (!argv$voom_off) {do_voom <- "--voom_off"} else {do_voom <- ""}
  if (!argv$ebayes_off) {do_ebayes <- "--ebayes_off"} else {do_ebayes <- ""}
  args <- paste("Rscript differential_gene_expression.r \\ \n  ",
    args$infile, " \\ \n  ", args$targets, " \\ \n  -o", args$outfile_dir,
    " \\ \n  -n", args$names, " \\ \n  ", do_voom, " \\ \n  ", do_ebayes, "\n"
  )
  print(cat(args))
  write(args, sep="", file=argpath)
}

main <- function() {
  argv <- parse_argv()
  dir.create(file.path(argv$outfile_dir))
  write_args(argv,paste(argv$outfile_dir,"/",basename(argv$names),".r",sep=""))
  pdf(paste(argv$outfile_dir, "/", basename(argv$names), ".pdf", sep=""))
  targets <- limma::readTargets(argv$targets)

  # create a design matrix, first column must be treatment condition
  celltype <- factor(targets[[1]])
  design <- model.matrix(~targets[[1]])

  # count numbers of reads mapped to NCBI Refseq genes
  data <- read.table(argv$infile, sep="\t")

  # for now any NA values are assumed to have 0 coverage
  data[is.na(data)] = 0
  dge <- edgeR::DGEList(counts=data)

  # filter out low-count genes
  isexpr <- rowSums(cpm(dge) > 10) >= 2
  dge <- dge[isexpr,]

  if (!argv$voom_off) {
    # perform voom normalization
    voomed <- limma::voom(dge, design, plot=TRUE)
  } else {
    voomed <- dge
  }

  # cluster libraries
  limma::plotMDS(voomed, xlim=c(-2.5,2.5))

  if (!argv$ebayes_off) {
    # fit linear model and assess differential expression
    fit <- limma::eBayes(limma::lmFit(voomed, design))
  } else {
    fit <- limma::lmFit(voomed, design)
  }

  write.table(limma::topTable(fit, number=Inf, sort.by="P", coef=2), sep="\t",
    file=paste(argv$outfile_dir, "/", basename(argv$names), ".top.tsv", sep=""),
    quote=F
  )
  dev.off()
}

main()
