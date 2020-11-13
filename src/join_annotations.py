#!/usr/bin/python
# take a list of counts file(s), join all into single file, rescaling if needed.
import argparse
import os
import subprocess as sp
from time import time
import warnings
import pandas as pd

def join_data(infile_path: str, annot_path: str, outfile_path: str=None):
    """
    Take a list of tab separated files with identical indices, join all.

    Arguments:
        (REQUIRED) infile_path: input data file path
        (REQUIRED) annot_path: input annotations file path
        (OPTIONAL) outfile_path: output file path

    Example input:

        input:
                   SAMPLE_1 ...
        FEATURE_1      1123 ...
        FEATURE_2       110 ...

        annotation:
                   SAMPLE_1_ANNOT   ...
        FEATURE_1             foo   ...
        FEATURE_2             bar   ...
        FEATURE_3            spam   ...

        output:
                   SAMPLE_1     SAMPLE_1_ANNOT
        FEATURE_1      1123                foo
        FEATURE_2       110                bar
    """
    data = pd.read_csv(infile_path, sep="\t", index_col=0)
    annot = pd.read_csv(annot_path, sep="\t", index_col=0)
    combined = pd.concat([data, annot], join="inner", axis=1)
    if outfile_path:
        combined.to_csv(outfile_path, sep="\t")
    return combined

def _argument_parser():
    parser = argparse.ArgumentParser(description=
        "Join tsv files with scores to annotations. \
        Usage: python join_annotations.py </path/to/data/*.tsv> \
            </path/to/annotations/*.tsv> -o [/path/to/out.tsv]"
        )
    parser.add_argument("infile_path", type=str,
                        help="Provide path to abundance file.")
    parser.add_argument("annot_path", type=str,
                        help="Provide path to annotation file.")
    parser.add_argument("-o", "--outfile_path", type=str, default=None,
                        help="Provide path to output file.")
    return parser.parse_args()

def main():
    args = _argument_parser()
    infile_path = args.infile_path
    annot_path = args.annot_path
    outfile_path = args.outfile_path

    print("# Loading data:")
    print(infile_path)
    print("# Loading annotations:")
    print(annot_path)

    if outfile_path is None:
        outfile_path = ".".join(["reannotated", str(int(time())), "tsv"])
        print("# No output file specified, writing to input dir:", outfile_path)

    if os.path.exists(outfile_path):
        warnings.warn("# Output file exists, overwriting!",category=UserWarning)
        os.remove(outfile_path)

    data = join_data(infile_path, annot_path, outfile_path)

    sp.call("".join(["perl -pi -e \'s/^\t//\' ", outfile_path]), shell=True)
    print("# Writing file to:", outfile_path)
    print("# Reproduce by running this command:")
    print("python join_annotations.py",infile_path,annot_path,"-o",outfile_path)

if __name__ == "__main__":
    main()
