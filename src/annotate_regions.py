#!/usr/bin/python
# take a matching gtf to genome regions file and adds annotations
import argparse
import os
import re
import subprocess as sp
from time import time
import warnings
import pandas as pd
from tqdm import tqdm

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
    return(regions[regions[filter_col] <= filter_val])

def annotate_regions(abundance: pd.DataFrame, gtf: pd.DataFrame,
                     outfile_path: str="./out.tsv", hide_progress: bool=True):
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
    # TODO: need to optimise, change strategy of matching overlaps?
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
                        help="Provide path to output file.")
    parser.add_argument("-p", "--hide_progress", action="store_true",
                        help="Hide the progress bar.")
    return parser.parse_args()

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

    annotate_regions(filter_regions(load_regions(infile_path),"adj.P.Val",0.05),
                     load_gtf(gtf_path),
                     outfile_path,
                     hide_progress=args.hide_progress)

if __name__ == "__main__":
    main()
