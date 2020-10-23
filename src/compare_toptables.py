#!/usr/bin/python
# extract gtf-like annotations and intersect gene names
import argparse
import os
import re
import subprocess as sp
from time import time
import warnings
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
from upsetplot import from_contents, UpSet

def load_filter_regions(infile_path: str, filter_col: str="adj.P.Val",
                        filter_val: float=0.05):
    """
    Load and filter regions of interest only, default =< 0.05 adjusted p vals

    Arguments:
        (REQUIRED) regions: dataframe containing data of interest
        (OPTIONAL) filter_col: column to filter on
        (OPTIONAL) filter_val: value in column to filter on
    """
    regions = pd.read_csv(infile_path, sep="\t")
    return regions[regions[filter_col] <= filter_val]

def extract_genes(data: pd.DataFrame, colname: str=None) -> list:
    """
    Take dataframe as input, extract feature names

    Arguments:
        (REQUIRED) data: path to input dataframe
        (OPTIONAL) colname: column name, if absent uses index
    """
    if not colname:
        return pd.Series(data.index).unique()
    else:
        return pd.Series(data[colname]).unique()

def parse_gtf(data: pd.Series, extract_col: str="gene_id"):
    """
    Take the attribute column of the gtf file as input, parse out feature names.

    More information on the gtf file format is here:
        https://asia.ensembl.org/info/website/upload/gff.html

    Arguments:
        (REQUIRED) data: pandas single data series
        (OPTIONAL) extract_col: column to extract (None returns whole dataframe)
    """
    # print(data.iloc[0])
    data = data.str.split(";")
    data = [[y.strip() for y in x] for x in data]
    data = [[y.replace("\"", "") for y in x] for x in data]
    data = [[y for y in x if y] for x in data]
    data = [[y.split(" ", 1) for y in x] for x in data]
    data = [x for x in _parse_gtf(data)]
    data = pd.DataFrame(data)
    if extract_col:
        return data[extract_col].unique()
    else:
        return data

def _parse_gtf(data):
    """Generator for parsing gtf file line by line"""
    for entry in data:
        row = dict()
        for features in entry:
            row[features[0]] = features[1]
        yield pd.Series(row)

def plot_intersection(data: dict, plot_outfile: str="upsetplot.pdf"):
    """
    Take a dict of lists of unique identifiers, make quantitative venn diagram.

    Arguments:
        (REQUIRED) data: dict of lists, transformed with from_contents
        (OPTIONAL) plot_outfile: save the figure here
    """
    data = UpSet(data,
                 show_counts=True,
                 show_percentages=True,
                 sort_categories_by=None)
    data.plot()
    if plot_outfile:
        plt.savefig(plot_outfile)
    else:
        plt.show()

def get_intersection(data: dict, sets_outfile: str="upsetplot.tsv"):
    """
    Take a dict of lists of unique identifiers, make quantitative venn diagram.

    Arguments:
        (REQUIRED) data: dict of lists, values within each list must be unique.
        (OPTIONAL) sets_outfile: save the set memberships here
    """
    data = from_contents(data)
    if sets_outfile:
        memberships = data.reset_index().set_index("id")
        memberships["present_in_sets"] = memberships.sum(axis=1)
        memberships.to_csv(sets_outfile, sep="\t")
    return data

def _argument_parser():
    timestamp = str(int(time()))
    parser = argparse.ArgumentParser(description=
        """
        Compare and contrast the attributes of multiple differentially
        expressed gene lists. No input validation is performed!
        Usage: python compare_toptables.py </path/to/original.tsv>
            </path/to/data.tsv> ... -o [/path/to/out.pdf]
        """
    )
    parser.add_argument("reference_path", type=str,
                        help="Provide path to reference toptables file. \
                        Must have gene name or identifier as column or index!")
    parser.add_argument("infile_paths", type=str, nargs="+",
                        help="Provide path to other toptables files. \
                        Must have gtf-like attribute annotation fields!")
    parser.add_argument("-p", "--plot_outfile", type=str,
                        default=".".join(["upset", timestamp, "pdf"]),
                        help="Provide path to output image file [eps/pdf/png].")
    parser.add_argument("-s", "--sets_outfile", type=str,
                        default=".".join(["upset", timestamp, "tsv"]),
                        help="Provide path to output set membership files.")
    parser.add_argument("-t", "--threshold", type=float, default=0.05,
                        help="Filter on adjusted pval threshold.")
    return parser.parse_args()

def main():
    args = _argument_parser()
    ref_path = args.reference_path
    infile_paths = args.infile_paths
    plot_outfile = args.plot_outfile
    sets_outfile = args.sets_outfile

    print("# No input validation is performed!")
    print("# Processing files:")
    print("#", ref_path)
    [print("#", infile_path) for infile_path in infile_paths]

    if sets_outfile is None:
        print("# No output file specified, not writing output!")

    if os.path.exists(sets_outfile):
        warnings.warn("Output file exists, overwriting!", category=UserWarning)

    if plot_outfile is None:
        print("# No output plot specified, not writing output!")

    if os.path.exists(plot_outfile):
        warnings.warn("Output plot exists, overwriting!", category=UserWarning)

    print("# Reproduce by running this command:")
    print(" ".join(["python compare_toptables.py", ref_path,
                    " ".join(infile_paths), "-p", plot_outfile,
                    "-s", sets_outfile, "-t", str(args.threshold)]))

    ref_key = ref_path.split("/")[-1].split(".")[0]
    ref = load_filter_regions(ref_path, filter_val=args.threshold)
    ref = extract_genes(ref)
    data_keys = [i.split("/")[-1] for i in infile_paths]
    data_keys = [i.split(".")[0] for i in data_keys]
    data = [load_filter_regions(i, filter_val=args.threshold)["attribute"]
            for i in infile_paths]
    data = [parse_gtf(i) for i in data]

    to_intersect = dict(zip(data_keys, data))
    to_intersect[ref_key] = ref

    intersection = get_intersection(to_intersect, sets_outfile)
    plot_intersection(intersection, plot_outfile)

if __name__ == "__main__":
    main()
