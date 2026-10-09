"""Count how reads spread over fixed genome bins, the data behind a fingerprint.

Sums the + and - strand 5' count tracks into bins over every chromosome in
the track header, which is the allow-list. Writes how many bins hold each read
count. Sorting the bins by count and accumulating gives the fingerprint curve.
"""

import argparse

import numpy as np
import polars as pl
import pyBigWig


def build_parser():
    """Build the fingerprint.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the fingerprint.py flags.
    """
    p = argparse.ArgumentParser(description="Count reads per genome bin.")
    p.add_argument("--plus", required=True, help="+ strand 5' count bigWig.")
    p.add_argument("--minus", required=True, help="- strand 5' count bigWig.")
    p.add_argument("--bin-size", type=int, default=500, help="Bin size in bp.")
    p.add_argument("-o", "--output", required=True, help="TSV to write.")
    return p


def main(argv=None):
    """Bin both strands and write the bin count distribution.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    plus, minus = pyBigWig.open(args.plus), pyBigWig.open(args.minus)
    per_chrom = []
    for chrom, length in plus.chroms().items():
        counts = np.zeros(-(-length // args.bin_size), dtype=np.int64)
        for bw in (plus, minus):
            intervals = bw.intervals(chrom) or []
            intervals = np.array(intervals, dtype=np.int64).reshape(-1, 3)
            np.add.at(counts, intervals[:, 0] // args.bin_size, intervals[:, 2])
        per_chrom.append(counts)
    reads, bins = np.unique(np.concatenate(per_chrom), return_counts=True)
    pl.DataFrame({"reads": reads, "bins": bins}).write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()
