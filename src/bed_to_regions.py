#!/usr/bin/python
# take a matching gtf to bedfile and adds annotations
import argparse
import os
import re
import subprocess as sp
from time import time
import warnings
from intervaltree import Interval, IntervalTree
import pandas as pd
from tqdm import tqdm

def load_bed(infile_path: str):
    """
    Load and parse a bed file. More information on the bedfile format is here:
        https://genome.ucsc.edu/FAQ/FAQformat.html#format1

    Arguments:
        (REQUIRED) infile_path: path to input bed file
    """
    return pd.read_csv(infile_path, sep="\t", header=None, index_col=0)

def load_sizes(infile_path: str, header: bool=None):
    """
    Load and parse a gtf file. More information on the gtf format is here:
        https://asia.ensembl.org/info/website/upload/gff.html

    Arguments:
        (REQUIRED) infile_path: path to gtf file
        (OPTIONAL) header: headers in size file (DEFAULT: None)

    chr1	247249719
    chr2	242951149
    ...
    """
    return pd.read_csv(infile_path, sep="\t", header=None, index_col=0)

def bed_to_windows(bed: pd.DataFrame, outfile_path: str="./out.bed",
                   window_size: int=1000, hide_progress: bool=True,
                   columns: list=["chrom", "chromStart", "chromEnd",
                                  "name", "score", "chromSize"],
                   aggregate: str=None):
    """
    Associate a coordinate-based genomic feature to a window instead.

    Arguments:
        (REQUIRED) bed: abundance measurements in dataframe
        (OPTIONAL) outfile_path: path to output file (DEFAULT: ./out.tsv)
        (OPTIONAL) window_size: size of windows to slice genome into
        (OPTIONAL) hide_progress: show the progress bar (DEFAULT: True)
        (OPTIONAL) columns: columns should be a normal bed file + chromSize
        (OPTIONAL) aggregate: aggregate the score columns (DEFAULT: sum)
            same as pandas aggregate functions, eg: [first,sum,median,mean,None]

    Example input:

        chrom  chromStart  chromEnd                   name     score  chromSize
        chr1       703964    703965   ICRH79alpha2_peak_1a  17.06508  247249719
        chr1       704109    704110   ICRH79alpha2_peak_1b  12.58500  247249719

    Example output:

                   SAMPLE_1     SAMPLE_2   ...  Gene id  ...
        FEATURE_1      1123         3402   ...     FOO1  ...
        FEATURE_1      1123         3402   ...     FOO2  ...
        FEATURE_2       110          149   ...      BAR  ...
    """
    # create indices corresponding for reverse mapping later on
    bed["index"] = bed[columns[0]] + "_" + bed[columns[1]].astype(str) + \
                   "_" + bed[columns[2]].astype(str)
    bed.set_index("index", inplace=True)

    # dont want cross-chromosome overlaps, build one intervaltree per chromosome
    for chr in tqdm(bed[columns[0]].unique(), total=bed[columns[0]].nunique(),
                    desc="Annotating genomic regions by chromosome",
                    disable=hide_progress):
        feature_coords = IntervalTree()

        # index what we want to map to
        window_max = bed[bed[columns[0]] == chr]
        window_max = window_max[columns[-1]].unique()
        assert len(window_max) == 1, "Chromosome can only have one size!"
        window_max = window_max[0]

        leftover = window_max % window_size
        window_count = int(window_max / window_size)

        start = 0
        end = window_size
        indices = list()

        while end < window_max:
            index = "_".join([chr, str(start), str(end)])
            indices.append(index)
            start += window_size
            end += window_size

        # the last window may not be an even divisor of window size
        last_index = "_".join([chr, str(start), str(window_max)])
        indices.append(last_index)

        windows = pd.DataFrame(indices)

        windows[["chrom", "windowStart", "windowEnd"]] = \
            windows[0].str.split("_", expand=True)
        windows[["windowStart"]] = windows[["windowStart"]].astype(int)
        windows[["windowEnd"]] = windows[["windowEnd"]].astype(int)
        windows.rename(columns={0: "index"}, inplace=True)
        windows.set_index("index", inplace=True)

        mapped_features = dict()
        bed_chr = bed[bed[columns[0]] == chr]
        for bed_entry in bed_chr.iterrows():
            index = bed_entry[0]
            start = bed_entry[1][columns[1]]
            end = bed_entry[1][columns[2]]
            name = bed_entry[1][columns[3]]
            score = bed_entry[1][columns[4]]
            feature_coords[start:end] = (index, start, end, name, score)

        # map the coordinates to reference
        mapped_features = dict()
        for window in windows.iterrows():
            start = window[1]["windowStart"]
            end = window[1]["windowEnd"]
            index = window[0]
            features = feature_coords[start:end]
            # drop any which map to nothing
            if features:
                mapped_features[index] = [i[-1] for i in list(features)]

        # use the indices of data and reference to annotate data with references
        for i, j in mapped_features.items():
            key = pd.DataFrame(windows.loc[i]).T
            val = pd.DataFrame(j)
            # line based format, multiple annotations can map to the same region
            region = pd.concat([key]*len(val))
            annotated = pd.concat([region.reset_index(), val], axis=1)
            annotated.set_index("index", inplace=True)
            annotated.columns = ["chrom", "windowStart", "windowEnd",
                                 "featureId", "featureStart", "featureEnd",
                                 "featureName", "score"]
            # save memory by writing to disk per block, also checkpoints
            annotated.to_csv(outfile_path, sep="\t", mode="a", header=None)

    # add the headers back in
    data = pd.read_csv(
        outfile_path, header=None, sep="\t", low_memory=False
        ).drop(0, axis=1)
    data.columns = annotated.columns

    # sort by chromosome
    try:
        data_tmp = data.chrom.str.split("chr", expand=True)
        data_tmp[1] = data_tmp[1].astype(int)
        data_tmp.sort_values(1, inplace=True)
        data = data.reindex(data_tmp.index)
    except:
        warnings.warn("Chromosomes not sorted numerically",category=UserWarning)

    data["window"] = data["chrom"].astype(str) + "_" + \
        data["windowStart"].astype(str) + "_" + \
        data["windowEnd"].astype(str)
    data.set_index("window", inplace=True)

    if aggregate is not None:
        data = data.groupby(level=0).agg(
            {"chrom": "first", "windowStart": "first", "windowEnd": "first",
             "featureId": lambda x: ";".join(x.astype(str)),
             "featureStart": lambda x: ";".join(x.astype(str)),
             "featureEnd": lambda x: ";".join(x.astype(str)),
             "featureName": lambda x: ";".join(x.astype(str)),
             "score": "sum"}
            )
    data.to_csv(outfile_path, mode="w", sep="\t")#, index=None)
    return data

def merge_data_sizes(data: pd.DataFrame, sizes: pd.DataFrame,
                     columns: list=["chrom", "chromStart", "chromEnd",
                                    "name", "score", "chromSize"]):
    """Join two dataframes (bedfile, chromosome sizes) together (same index)."""
    data = pd.merge(data, sizes, left_index=True, right_index=True)
    data.reset_index(inplace=True)
    data.columns = columns
    return data

def _argument_parser():
    parser = argparse.ArgumentParser(description=
        "Transform a bed file to a coordinate window format. Can take gzip. \
        No input validation is performed! \
        Usage: python reformat_bedfile.py </path/to/data.tsv> </path/to/sizes> \
            -o [/path/to/out.tsv]"
        )
    parser.add_argument("infile_path", type=str,
                        help="Provide path to abundance measurements file.")
    parser.add_argument("chrom_sizes", type=str,
                        help="Provide path to chromosome sizes file.")
    parser.add_argument("-a", "--aggregate", type=str, default=None,
                        help="Aggregate scores [None,sum,median,mean,first]")
    parser.add_argument("-g", "--gtf_path", type=str, default=None,
                        help="Provide path to genome annotations file.")
    parser.add_argument("-o", "--outfile_path", type=str,
                        help="Provide path to output file (can be .gz).")
    parser.add_argument("-w", "--window_size", type=int, default=1000,
                        help="Window size to break genome into (DEFAULT: 1000)")
    parser.add_argument("-p", "--hide_progress", action="store_true",
                        help="Hide the progress bar.")
    return parser.parse_args()

def main():
    args = _argument_parser()
    infile_path = args.infile_path
    chrom_sizes = args.chrom_sizes
    outfile_path = args.outfile_path
    window_size = args.window_size

    print("# No input validation is performed!")
    print("# Make sure genome size file matches genome version.")
    print("# Processing file:")
    print("#", infile_path)
    print("# Using sizes:")
    print("#", chrom_sizes)

    if outfile_path is None:
        outfile_path = ".".join(["annotated", str(int(time())), "tsv"])
        print("# No output file specified, writing to current dir:",outfile_path)

    if os.path.exists(outfile_path):
        warnings.warn("# Output file exists, overwriting!",category=UserWarning)
        os.remove(outfile_path)

    print("# Reproduce by running this command:")
    if args.aggregate is not None:
        print(" ".join(["python bed_to_regions.py", infile_path, chrom_sizes,
                        "-a", args.aggregate, "-o", outfile_path,
                        "-w", str(window_size), "-p"]))
    else:
        print(" ".join(["python bed_to_regions.py", infile_path, chrom_sizes,
                        "-o", outfile_path, "-w", str(window_size), "-p"]))

    data = load_bed(infile_path)
    sizes = load_sizes(chrom_sizes)
    data_sizes = merge_data_sizes(data, sizes)
    windows = bed_to_windows(
        data_sizes, outfile_path=outfile_path, aggregate=args.aggregate,
        window_size=window_size, hide_progress=args.hide_progress
        )

if __name__ == "__main__":
    main()
