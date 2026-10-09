"""Average the + and - strand 5' counts around peak summits.

In ChIP-nexus, + strand 5' ends stop left of the bound protein and - strand
5' ends stop right of it. The averaged profile shows that footprint.
"""

import argparse

import numpy as np
import polars as pl
import pyBigWig
from zlib_ng import gzip_ng


def build_parser():
    """Build the footprint.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the footprint.py flags.
    """
    p = argparse.ArgumentParser(description="Average 5' counts around summits.")
    p.add_argument("--peaks", required=True, help="Gzipped narrowPeak with summits.")
    p.add_argument("--samples", nargs="+", required=True, help="Sample names.")
    p.add_argument(
        "--plus", nargs="+", required=True, help="+ strand bigWig per sample."
    )
    p.add_argument(
        "--minus", nargs="+", required=True, help="- strand bigWig per sample."
    )
    p.add_argument("--width", type=int, default=100, help="Half window in bp.")
    p.add_argument("-o", "--output", required=True, help="TSV to write.")
    return p


def main(argv=None):
    """Average each sample's tracks around the summits and write the profiles.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with gzip_ng.open(args.peaks, "rt") as f:
        summits = [
            (r[0], int(r[1]) + int(r[9])) for r in (line.split("\t") for line in f)
        ]

    offsets = np.arange(-args.width, args.width + 1)
    tables = []
    for sample, plus, minus in zip(args.samples, args.plus, args.minus):
        profiles = [_mean_profile(path, summits, args.width) for path in (plus, minus)]
        tables.append(
            pl.DataFrame(
                {
                    "sample": sample,
                    "offset": offsets,
                    "plus": profiles[0],
                    "minus": profiles[1],
                }
            )
        )
    pl.concat(tables).write_csv(args.output, separator="\t")


def _mean_profile(path, summits, width):
    """Average one track over windows centered on the summits.

    Parameters
    ----------
    path : str
        The bigWig of 5' counts.
    summits : list of tuple
        `(chrom, position)` per summit.
    width : int
        Half window in bp.

    Returns
    -------
    numpy.ndarray
        Mean count per offset. Zeros without usable summits.
    """
    bw = pyBigWig.open(path)
    sizes = bw.chroms()
    total, n = np.zeros(2 * width + 1), 0
    for chrom, pos in summits:
        if chrom in sizes and width <= pos < sizes[chrom] - width:
            total += np.nan_to_num(
                bw.values(chrom, pos - width, pos + width + 1, numpy=True)
            )
            n += 1
    bw.close()
    return total / n if n else total


if __name__ == "__main__":
    main()
