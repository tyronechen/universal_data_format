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
                      in column InputFile and OutputFile! If paired-end reads, \
                      include column InputFile2.")
    p <- add_argument(p, "base_name", type="character",
                      help="Base name for index, eg human_hg19.")
    p <- add_argument(p, "reference", type="character",
                      help="Reference genome as a fasta file [fa/fna/fasta]. \
                      This creates the index in that directory.")
    p <- add_argument(p, "--gtf_annot", type="character",
                      default=NULL,
                      help="Custom annotations file in gtf or gtf.gz format.")
    p <- add_argument(p, "--featuretype", type="character",
                      default="exon",
                      help="Annotate this feature from gtf.")
    p <- add_argument(p, "--featurecount_multimap", flag=TRUE,
                      help="Include multi-mapping reads in feature counting.")
    p <- add_argument(p, "--counts_path", type="character",
                      default="./counts.tsv",
                      help="Output file for abundance measurements.")
    p <- add_argument(p, "--assembly", type="character", default="hg19",
                      help="Genome assembly version [mm9, mm10, hg19, hg38].")
    p <- add_argument(p, "--nthreads", type="integer", default=1,
                      help="Number of threads to use.")
    p <- add_argument(p, "--mem", type="integer", default=16000,
                      help="Amount of memory (MB) to use.")
    argv <- parse_args(p)
    return(argv)
}

write_args <- function(args, argpath) {
  if ( args$featurecount_multimap == TRUE ) {
    args <- paste("Rscript align.r \\ \n  ",
      args$targets, " \\ \n  ", args$base_name, " \\ \n  ",
      args$reference, " \\ \n  -g", args$gtf_annot, " \\ \n  -c",
      args$counts_path, " \\ \n  -a", args$assembly, "\\ \n -t",
      args$featuretype, " \\ \n  -n", args$nthreads, " \\ \n  -m",
      args$mem, " \\ \n  --featurecount_multimap"
    )
  } else {
    args <- paste("Rscript align.r \\ \n  ",
      args$targets, " \\ \n  ", args$base_name, " \\ \n  ",
      args$reference, " \\ \n  -g", args$gtf_annot, " \\ \n  -c",
      args$counts_path, " \\ \n  -a", args$assembly, "\\ \n -t",
      args$featuretype, " \\ \n  -n", args$nthreads, " \\ \n  -m",
      args$mem, " \n"
    )
  }
  print(cat(args))
  write(args, sep="", file=argpath)
}

main <- function() {
  argv <- parse_argv()
  dir.create(file.path(dirname(argv$counts_path)))
  write_args(argv, paste(argv$counts_path, ".r", sep=""))

  print("foo")

  # read in target file
  options(digits=2)
  targets <- readTargets(argv$targets)

  # build an index for reference sequence (Chr1 in hg19)
  original_dir <- getwd()
  setwd(dirname(argv$reference))

  # if index exists
  suffixes <- paste(
    argv$base_name, c("00.b.array", "00.b.tab", "reads", "files"), sep="."
  )
  if ( !all(file.exists(suffixes)) ) {
    print("# No subread index matching that name found, building new index...")
    buildindex(
      basename=argv$base_name,
      reference=basename(argv$reference),
      memory=argv$mem
    )
  }

  # if reads exist, dont need to realign
  if ( !all(file.exists(targets$OutputFile)) ) {
    # align reads
    align(
      index=argv$base_name,
      readfile1=targets$InputFile,
      readfile2=targets$InputFile2,
      input_format="gzFASTQ",
      output_format="BAM",
      output_file=targets$OutputFile,
      unique=TRUE,
      indels=5,
      nthreads=argv$nthreads,
    )
  } else {
    print("# Found files matching names in OutputFile column of targets file.")
    print("# Skipping alignment. Move or rename these files for new alignment.")
  }

  # paired end reads are quantified differently in featureCounts
  if ( !is.null(targets$InputFile) & !is.null(targets$InputFile2) ) {
    isPairedEnd <- TRUE
  } else {isPairedEnd <- FALSE}

  # inbuilt subread annotations only include a few versions
  if ( !argv$assembly %in% c("mm9", "mm10", "hg19", "hg38") ) {
    argv$assembly <- NULL
  }

  # count numbers of reads mapped to NCBI Refseq genes
  fc <- featureCounts(
    files=targets$OutputFile,
    annot.inbuilt=argv$assembly,
    annot.ext=argv$gtf_annot,
    isGTFAnnotationFile=TRUE,
    GTF.featureType=argv$featuretype,
    isPairedEnd=isPairedEnd,
    countMultiMappingReads=argv$featurecount_multimap,
    nthreads=argv$nthreads
  )
  setwd(original_dir)
  write.table(fc$counts, sep="\t", quote=F,
    file=paste(argv$outfile_dir, "/", argv$counts_path, ".tsv", sep="")
  )
}

main()
