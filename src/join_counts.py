#!/usr/bin/python
# take a list of counts file(s), join all into single file, rescaling if needed.
import argparse
import os
import re
import subprocess as sp
from time import time
import warnings
import numpy as np
import pandas as pd

def join_data(infile_paths: list, outfile_path: str, rescale: bool=False):
    """
    Take a list of tab separated files with identical indices, join all.

    Arguments:
        (REQUIRED) infile_paths: list of input file paths
        (REQUIRED) outfile_path: output file path
        (OPTIONAL) rescale: if True, rescale all values to a scale of [0, 1]

    Example input:

        file 1:
                   SAMPLE_1
        FEATURE_1      1123
        FEATURE_2       110

        file 2:
                   SAMPLE_2
        FEATURE_1      3402
        FEATURE_2       149

        output:
                   SAMPLE_1     SAMPLE_2
        FEATURE_1      1123         3402
        FEATURE_2       110          149
    """
    all_samples = [pd.read_csv(i, sep="\t", index_col=0) for i in infile_paths]
    all_samples = pd.concat(all_samples, axis=1)
    if rescale:
        all_samples -= all_samples.min()
        all_samples /= all_samples.max()
        is_empty = all_samples[(all_samples < 0) | (all_samples > 1)].dropna()
        assert is_empty.empty, "All values should be [0,1]! Check for negative?"
    return all_samples

def join_contiguous(data: pd.DataFrame, filter_val: int=0):
    """
    Join adjacent genomic regions together where reads map.

    Arguments:
        (REQUIRED) data: genomic regions labelled as <chrname>_<start>_<end>
        (OPTIONAL) filter_val: dont join regions where value is sparse

    Note that entries across all samples will be kept if one element is detected
    You may want to filter small values before this step if they are not needed

    Example input:
        input:
                        SAMPLE_1     SAMPLE_2   ...
        chr1_001_100           0            0   ...
        chr1_101_200         100           90   ...
        chr1_201_300          50           70   ...
        chr1_301_400           0           20   ...
        chr1_401_500           0            0   ...
        chr1_501_600           0            0   ...

        output:
                        SAMPLE_1     SAMPLE_2   ...
        chr1_001_100           0            0   ...
        chr1_101_400         150          180   ...
        chr1_401_600           0            0   ...

    """
    is_contig = ((data > filter_val).sum(axis=1) > 0)
    not_contig = ((data <= filter_val).sum(axis=1) == data.shape[1])

    contigs = np.r_[False, is_contig, False]
    contigs = np.diff(contigs).nonzero()[0]
    contigs = np.reshape(contigs, (-1,2))

    gaps = np.r_[False, not_contig, False]
    gaps = np.diff(gaps).nonzero()[0]
    gaps = np.reshape(gaps, (-1,2))

    contigs = [combine_contigs(data, i) for i in contigs]
    gaps = [combine_contigs(data, i) for i in gaps]

    data = pd.concat([pd.concat(contigs), pd.concat(gaps)], axis=0)
    new_indices = [i[0] if type(i) == tuple else i for i in data.index]
    data.index = new_indices
    return data.sort_index()

def combine_contigs(data: pd.DataFrame, contig_indices: np.ndarray):
    """
    Add elements in sample column to match contiguous region

    Arguments:
        (REQUIRED) data: pandas dataframe to be collapsed
        (REQUIRED) contig_indices: indices of contigs to take
    """
    data = data.iloc[contig_indices[0]:contig_indices[1]]
    index = data.index
    if len(index) > 1:
        start = index[0].split("_")
        end = index[-1].split("_")
        index = "_".join([start[0], start[1], end[2]])
    data = data.sum(axis=0)
    data = pd.DataFrame(data).T
    data.index = [index]
    return data

def _argument_parser():
    parser = argparse.ArgumentParser(description=
        "Join tsv files corresponding to individual samples. \
        Usage: python join_counts.py </path/to/data/*.tsv> ... \
            -o [/path/to/out.tsv]"
        )
    parser.add_argument("infile_paths", type=str, nargs='+',
                        help="Provide path to multiple abundance files.")
    parser.add_argument("-o", "--outfile_path", type=str,
                        help="Provide path to output file.")
    parser.add_argument("-r", "--rescale", action="store_true",
                        help="Rescale values to a scale of [0,1].")
    parser.add_argument("-m", "--make_contiguous", type=int, default=None,
                        help="Join adjacent genomic regions <= filter.")
    return parser.parse_args()

def main():
    args = _argument_parser()
    infile_paths = args.infile_paths
    outfile_path = args.outfile_path

    print("# Processing files:")
    [print("#", infile) for infile in infile_paths]

    if outfile_path is None:
        outfile_path = ".".join(["all_samples", str(int(time())), "tsv"])
        print("# No output file specified, writing to input dir:", outfile_path)

    if os.path.exists(outfile_path):
        warnings.warn("# Output file exists, overwriting!", category=UserWarning)
        os.remove(outfile_path)

    data = join_data(infile_paths, outfile_path, rescale=args.rescale)
    if not args.make_contiguous is None:
        data = join_contiguous(data, filter_val=args.make_contiguous)

    data.to_csv(outfile_path, sep="\t")
    sp.call("".join(["perl -pi -e \'s/^\t//\' ", outfile_path]), shell=True)
    print("# Writing file to:", outfile_path)

    print("# Reproduce by running this command:")
    if args.rescale:
        print("python join_counts.py", " ".join(infile_paths),
              "-o", outfile_path, "-r", "-m", args.make_contiguous)
    else:
        print("python join_counts.py", " ".join(infile_paths),
              "-o", outfile_path, "-m", args.make_contiguous)

if __name__ == "__main__":
    main()
