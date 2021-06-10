#!/usr/bin/python
# take a list of counts file(s), join all into single file, rescaling if needed.
import argparse
import json
import os
import re
import subprocess as sp
from time import time
import warnings
# from numba import njit, prange, set_num_threads
import numpy as np
import pandas as pd

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
                        help="Join adjacent genomic regions <= filter.\
                        (not yet supported, will be added in future)")
    # parser.add_argument("-n", "--ncpus", type=int, default=1,
    #                     help="Specify number of cpus for parallelising \
    #                     operations. (default 1)")
    parser.add_argument("-s", "--split", type=str, default=None,
                        help="Split scores from annotations. Useful if data is \
                        already annotated and you want to keep the annotation. \
                        If this option is enabled, annotations will be written \
                        to outfile_path.annot. Sample argument: \
                        \'{\"names\": [\"sample1\", ... ], \
                        \"keep\":[\"col1\", ... ]}\'")
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

def join_data(infile_paths: list, outfile_path: str, rescale: bool=False,
              columns: dict=None):
    """
    Take a list of tab separated files with identical indices, join all.

    Arguments:
        (REQUIRED) infile_paths: list of input file paths
        (REQUIRED) outfile_path: output file path
        (OPTIONAL) rescale: if True, rescale all values to a scale of [0, 1]
        (OPTIONAL) columns: keep these columns.
            Example: {keep: ["col1", "col2", ... ]
                      names: ["sample1", "sample2", ... ]}

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
    data = [pd.read_csv(i, sep="\t", index_col=0) for i in infile_paths]
    data = pd.concat(data, axis=1)

    if columns:
        if "names" in columns:
            names = columns["names"]
            names = [str(x) for x in names]
        else:
            names = [os.path.basename(i) for i in infile_paths]
            names = [i.split("_")[0] for i in names if i.startswith("GSM")]
        keep = columns["keep"]
        assert len(infile_paths) == len(names), \
            "Number of sample names must match number of input files!"
        data = data[keep]

        # the first column series contains scores, set sample name as header
        newcols = [["_".join([y, x]) for y in names] for x in keep[1:]]
        newcols = [x for y in newcols for x in y]
        newcols = names + newcols
        data.columns = newcols

        # counts and annotations should be separate, row ids same for rejoining
        outfile_path_annot = ".".join([outfile_path, "annot"])
        annots = pd.concat([data.filter(regex=i, axis=1) for i in keep[1:]])
        annots.to_csv(outfile_path_annot, sep="\t")

        # will pass to R which requires a highly specific dataframe input format
        sp.call("".join(["perl -pi -e \'s/^\t//\' ", outfile_path_annot]), shell=True)
        print("# Writing annotations file to:", outfile_path_annot)
        data = data[names]

    if rescale:
        data -= data.min()
        data /= data.max()
        is_empty = data[(data < 0) | (data > 1)].dropna()
        assert is_empty.empty, "All values should be [0,1]! Check for negative?"
    return data

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
    # TODO: placeholder, nonfunctional!
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
    data["sort_key"] = data.index #.str.split("_", expand=True)
    data[["seqnames", "start", "end"]] = data.sort_key.str.split("_", expand=True)
    data["start"] = data["start"].astype(int)
    data.sort_values("start", inplace=True)
    data.drop(["sort_key", "seqnames", "start", "end"], axis=1, inplace=True)
    return data

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

def main():
    args = _argument_parser()
    infile_paths = args.infile_paths
    outfile_path = args.outfile_path
    columns = args.split

    print("# Processing files:")
    [print("#", infile) for infile in infile_paths]

    if outfile_path is None:
        outfile_path = ".".join(["all_samples", str(int(time())), "tsv"])
        print("# No output file specified, writing to input dir:", outfile_path)

    if os.path.exists(outfile_path):
        warnings.warn("# Output file exists, overwriting!", category=UserWarning)
        os.remove(outfile_path)

    if columns:
        columns = json.loads(columns)
    data = join_data(infile_paths, outfile_path, rescale=args.rescale,
                     columns=columns)

    if not args.make_contiguous is None:
        data = join_contiguous(data, filter_val=args.make_contiguous)

    data.to_csv(outfile_path, sep="\t")
    sp.call("".join(["perl -pi -e \'s/^\t//\' ", outfile_path]), shell=True)
    print("# Writing file to:", outfile_path)

    print("# Reproduce by running this command:")
    if args.rescale:
        print("python join_counts.py", " ".join(infile_paths),
              "-o", outfile_path, "-r", "-m", args.make_contiguous,
              "-s", columns)
    else:
        print("python join_counts.py", " ".join(infile_paths),
              "-o", outfile_path, "-m", args.make_contiguous,
              "-s", columns)

if __name__ == "__main__":
    main()
