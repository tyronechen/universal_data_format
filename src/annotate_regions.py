#!/usr/bin/python
# take a matching gtf to genome regions file and adds annotations
import argparse
import os
import re
import subprocess as sp
from time import time
import warnings
from intervaltree import Interval, IntervalTree
# from numba import njit, prange, set_num_threads
import pandas as pd
from tqdm import tqdm

def _argument_parser():
    parser = argparse.ArgumentParser(description=
        "Annotate a table of genomic regions with a gtf file. Can take gzip. \
        No input validation is performed! \
        Usage: python annotate_regions.py </path/to/data.tsv> </path/to/gtf> \
            -o [/path/to/out.tsv]"
        )
    parser.add_argument("infile_path", type=str,
                        help="Provide path to abundance measurements file.")
    parser.add_argument("gtf_path", type=str,
                        help="Provide path to genome annotations file.")
    parser.add_argument("-o", "--outfile_path", type=str,
                        help="Provide path to output file (can be .gz).")
    parser.add_argument("-t", "--threshold", type=float, default=0.05,
                        help="Filter on adjusted pval threshold (if DEG list).")
    # parser.add_argument("-n", "--ncpus", type=int, default=1,
    #                     help="Specify number of cpus for parallelising \
    #                     operations. (default 1)")
    parser.add_argument("-p", "--hide_progress", action="store_true",
                        help="Hide the progress bar.")
    parser.add_argument("-l", "--lowmem", action="store_true",
                        help="Very slow annotation but requires less memory. \
                        Dont use this unless necessary, O(m x n) from O(n)!")
    return parser.parse_args()

def check_parallel():
    args = _argument_parser()
    infile_paths = args.infile_paths
    if  args.ncpus > 1:
        do_parallel = True
    else:
        do_parallel = False
    return args.ncpus, do_parallel

# ncpus, do_parallel = check_parallel()
# set_num_threads(ncpus)

def load_regions(infile_path: str):
    """
    Load and parse a tab separated file of abundance measurements/dge. As long
    as the chromosome id, start indices and end indices are preserved, the
    content can be anything. Must have a header line.

    Arguments:
        (REQUIRED) infile_path: path to input file

    Example input:
        index           anything can go here    any number of columns
        chr1_1000_10000                  ...                      ...
    """
    data = pd.read_csv(infile_path, sep="\t", index_col=0)
    data.reset_index(inplace=True)
    data.rename(columns={"Unnamed: 0": "index"}, inplace=True)
    data.columns.str.replace(r'.bam$', '')
    data[["region_names", "start_region", "end_region"]] = \
        data["index"].str.split("_", expand=True,)
    data[["start_region"]] = data[["start_region"]].astype(int)
    data[["end_region"]] = data[["end_region"]].astype(int)
    return data

def load_gtf(infile_path: str, header: list=["seqnames", "source", "feature",
             "start", "end", "score", "strand", "frame", "attribute"]):
    """
    Load and parse a gtf file. More information on the gtf format is here:
        https://asia.ensembl.org/info/website/upload/gff.html

    Arguments:
        (REQUIRED) infile_path: path to gtf file
        (OPTIONAL) header: headers in gtf file (DEFAULT: standard gtf headers)
    """
    return pd.read_csv(infile_path, sep="\t", skiprows=5, names=header)

def filter_regions(regions: pd.DataFrame, filter_col: str="adj.P.Val",
                   filter_val: float=0.05):
    """
    Filter out regions of interest only, default =< 0.05 adjusted p vals

    Arguments:
        (REQUIRED) regions: dataframe containing data of interest
        (OPTIONAL) filter_col: column to filter on
        (OPTIONAL) filter_val: value in column to filter on
    """
    if filter_col in regions.columns:
        return regions[regions[filter_col] <= filter_val]
    else:
        warnings.warn("Column not found, not filtering!", category=UserWarning)
        return regions

def annotate_regions_lowmem(abundance: pd.DataFrame, gtf: pd.DataFrame,
                            outfile_path: str="./out.tsv",
                            hide_progress: bool=True):
    """
    Annotate region of genome with identifiers. The output format is line-based
    and multi-mapping is possible, for example if a region of interest contains
    more than one genomic feature. In these cases, the line will be duplicated!
    (No good reason to use this slow function, but leaving here in case needed)

    Arguments:
        (REQUIRED) abundance: abundance measurements in dataframe
        (REQUIRED) gtf: genome annotations in dataframe
        (OPTIONAL) outfile_path: path to output file (DEFAULT: ./out.tsv)
        (OPTIONAL) hide_progress: show the progress bar (DEFAULT: True)

    Example input (abundance):

                   SAMPLE_1     SAMPLE_2    ...
        FEATURE_1      1123         3402    ...
        FEATURE_2       110          149    ...

    Example input (gtf):

        seqnames    ...
        FOO         ...
        BAR         ...

    Example output:

                   SAMPLE_1     SAMPLE_2   ...  Gene id  ...
        FEATURE_1      1123         3402   ...     FOO1  ...
        FEATURE_1      1123         3402   ...     FOO2  ...
        FEATURE_2       110          149   ...      BAR  ...
    """
    for abu_entry in tqdm(abundance.iterrows(), total=abundance.shape[0],
                          desc="Annotating genomic regions",
                          disable=hide_progress):
        abu_start = abu_entry[1]["start_region"]
        abu_end = abu_entry[1]["end_region"]
        abu_chr = abu_entry[1]["region_names"]

        for gtf_entry in gtf.iterrows():
            gtf_start = gtf_entry[1]["start"]
            gtf_end = gtf_entry[1]["end"]
            gtf_chr = gtf_entry[1]["seqnames"]

            if abu_chr == gtf_chr:
                if len(range(max(gtf_start,abu_start),min(gtf_end,abu_end)+1))>0:
                    mapped = [abu_entry[1], gtf_entry[1]]
                    mapped = pd.DataFrame(pd.concat(mapped, axis=0)).transpose()
                    mapped.to_csv(outfile_path, sep="\t", mode="a", header=None)

    data = pd.read_csv(outfile_path, header=None, sep="\t").drop(0, axis=1)
    data.columns = mapped.columns
    data.to_csv(outfile_path, mode="w", sep="\t", index=None)
    return data

def annotate_regions_highmem(abundance: pd.DataFrame, gtf: pd.DataFrame,
                             outfile_path: str="./out.tsv",
                             hide_progress: bool=True, sort_by: str="adj.P.Val"):
    """
    Annotate region of genome with identifiers. The output format is line-based
    and multi-mapping is possible, for example if a region of interest contains
    more than one genomic feature. In these cases, the line will be duplicated!

    Arguments:
        (REQUIRED) abundance: abundance measurements in dataframe
        (REQUIRED) gtf: genome annotations in dataframe
        (OPTIONAL) outfile_path: path to output file (DEFAULT: ./out.tsv)
        (OPTIONAL) hide_progress: show the progress bar (DEFAULT: True)

    Example input (abundance):

                   SAMPLE_1     SAMPLE_2    ...
        FEATURE_1      1123         3402    ...
        FEATURE_2       110          149    ...

    Example input (gtf):

        seqnames    ...
        FOO         ...
        BAR         ...

    Example output:

                   SAMPLE_1     SAMPLE_2   ...  Gene id  ...
        FEATURE_1      1123         3402   ...     FOO1  ...
        FEATURE_1      1123         3402   ...     FOO2  ...
        FEATURE_2       110          149   ...      BAR  ...
    """
    # create indices corresponding for reverse mapping later on
    gtf["index"] = gtf["seqnames"] + "_" + gtf["start"].astype(str) + "_" + \
                   gtf["end"].astype(str)
    gtf.set_index("index", inplace=True)
    abundance.set_index("index", inplace=True)

    # intervaltree slicing is not inclusive of end point
    gtf["end"] = gtf["end"] + 1
    abundance["end_region"] = abundance["end_region"] + 1

    # dont want cross-chromosome overlaps, build one intervaltree per chromosome
    for chr in tqdm(gtf.seqnames.unique(), total=gtf.seqnames.nunique(),
                    desc="Annotating genomic regions by chromosome",
                    disable=hide_progress):
        ref_coords = IntervalTree()

        # index the reference chromosome
        gtf_chr = gtf[gtf["seqnames"] == chr]
        for gtf_entry in gtf_chr.iterrows():
            start = gtf_entry[1]["start"]
            end = gtf_entry[1]["end"]
            index = gtf_entry[0]
            ref_coords[start:end] = index

        # map the coordinates to reference
        mapped_features = dict()
        abu_chr = abundance[abundance["region_names"] == chr]
        for abu_entry in abu_chr.iterrows():
            start = abu_entry[1]["start_region"]
            end = abu_entry[1]["end_region"]
            index = abu_entry[0]
            features = ref_coords[start:end]
            # drop any which map to nothing
            if features:
                mapped_features[index] = [i[-1] for i in list(features)]

        # use the indices of data and reference to annotate data with references
        for i, j in mapped_features.items():
            key = pd.DataFrame(abundance.loc[i]).T
            val = gtf.loc[j].reset_index()
            # line based format, multiple annotations can map to the same region
            region = pd.concat([key]*len(val))
            val.rename(columns={"index": "index_reference"}, inplace=True)
            annotated = pd.concat([region.reset_index(), val], axis=1)
            annotated["end"] = annotated["end"] - 1
            annotated["end_region"] = annotated["end_region"] - 1
            # save memory by writing to disk per block, also checkpoints
            annotated.to_csv(outfile_path, sep="\t", mode="a", header=None)

    # add the headers back in
    data = pd.read_csv(outfile_path, header=None, sep="\t").drop(0, axis=1)
    data.columns = annotated.columns
    if sort_by in data.columns:
        data.sort_values(sort_by, ascending=True, inplace=True)
    else:
        print("# No sorting of values performed.")
    data.to_csv(outfile_path, mode="w", sep="\t", index=None)
    return data

def main():
    args = _argument_parser()
    infile_path = args.infile_path
    gtf_path = args.gtf_path
    outfile_path = args.outfile_path

    print("# No input validation is performed!")
    print("# Make sure gtf file is formatted correctly, matches genome version.")
    print("# Processing file:")
    print("#", infile_path)
    print("# Using reference:")
    print("#", gtf_path)

    if outfile_path is None:
        outfile_path = ".".join(["annotated", str(int(time())), "tsv"])
        print("# No output file specified, writing to current dir:",outfile_path)

    if os.path.exists(outfile_path):
        warnings.warn("# Output file exists, overwriting!",category=UserWarning)
        os.remove(outfile_path)

    print("# Reproduce by running this command:")
    print(" ".join(["python annotate_regions.py", infile_path, gtf_path,
                    "-o", outfile_path, "-t", str(args.threshold), "-p"]))

    data = load_regions(infile_path)

    if args.threshold:
        data = filter_regions(data, "adj.P.Val", args.threshold)

    if args.lowmem is False:
        annotate_regions_highmem(data,
                                 load_gtf(gtf_path),
                                 outfile_path,
                                 hide_progress=args.hide_progress)
    else:
        annotate_regions_lowmem(data,
                                load_gtf(gtf_path),
                                outfile_path,
                                hide_progress=args.hide_progress)

if __name__ == "__main__":
    main()
