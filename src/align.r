#!/usr/bin/Rscript
# This is a modified version of the code found here:
#   http://bioinf.wehi.edu.au/RNAseqCaseStudy/

# To run this case study, you should have R version of 3.0.2 or later

# load libraries
library(Rsubread)
library(limma)
library(edgeR)

parse_argv <- function() {
    library(argparser, quietly=TRUE)
    p <- arg_parser("Align reads into bam files and export abundance measures.")
    p <- add_argument(p, "targets", type="character",
                      help="Sample information. Must include file paths \
                      in column InputFile and OutputFile!")
    p <- add_argument(p, "basename", type="character",
                      help="Base name for index, eg human_hg19")
    p <- add_argument(p, "reference", type="character",
                      help="Reference genome as a fasta file [fa/fna/fasta]. \
                      This creates the index in that directory.")
    p <- add_argument(p, "--counts_path", type="character",
                      default="./counts.tsv",
                      help="Output file for abundance measurements")
    p <- add_argument(p, "--assembly", type="character", default="hg19",
                      help="Genome assembly version [hg19, hg38]")
    argv <- parse_args(p)
    return(argv)
}

write_args <- function(args, argpath) {
  args <- paste("Rscript differential_gene_expression.r \\ \n  ",
    args$targets, " \\ \n  ", args$basename, " \\ \n  ",
    args$reference, " \\ \n  -c", args$counts_path, " \\ \n  -a",
    args$assembly, "\n"
  )
  print(cat(args))
  write(args, sep="", file=argpath)
}

main <- function() {
  argv <- parse_argv()
  dir.create(file.path(dirname(argv$counts_path)))
  write_args(argv, paste(argv$counts_path, ".r", sep=""))

  # read in target file
  options(digits=2)
  targets <- readTargets(argv$targets)

  # build an index for reference sequence (Chr1 in hg19)
  original_dir <- getwd()
  setwd(dirname(argv$reference))
  buildindex(basename=argv$basename, reference=basename(argv$reference))
  setwd(original_dir)

  # align reads
  align(
    index=argv$basename,
    readfile1=targets$InputFile,
    input_format="gzFASTQ",
    output_format="BAM",
    output_file=targets$OutputFile,
    unique=TRUE,
    indels=5
  )

  # count numbers of reads mapped to NCBI Refseq genes
  fc <- featureCounts(files=targets$OutputFile, annot.inbuilt=argv$assembly)
  write.table(fc$counts, sep="\t", quote=F,
    file=paste(argv$outfile_dir, "/", argv$counts_path, ".tsv", sep="")
  )
}

main()
