#!/usr/bin/python
# take a mauve multiseq alignment, remove gaps & N and split to multifasta file.
import argparse
import os
import re
from time import time
import warnings
import pandas as pd

def join_data(infile_paths: list, outfile_path: str, rescale: bool=False):
    """Take a list of tab separated files with identical indices, join all.

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
    all_samples.to_csv(outfile_path, sep="\t")
    return all_samples

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

    join_data(infile_paths, outfile_path, rescale=args.rescale)
    
    print("# Writing file to:", outfile_path)

    print("# Reproduce by running this command:")
    if args.rescale:
        print("python join_counts.py", " ".join(infile_paths),
              "-o", outfile_path, "-r")
    else:
        print("python join_counts.py"," ".join(infile_paths),"-o",outfile_path)

if __name__ == "__main__":
    main()
