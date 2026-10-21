"""Plot the fingerprints of a ChIP experiment and its control.

Bins are sorted by read count. The curve is the fraction of reads in the
lowest-count fraction of bins. Uniform coverage follows the diagonal. A strong
ChIP stays low and rises late, because a few bins hold most reads.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import poikilos as pk
import polars as pl


def build_parser():
    """Build the plot_fingerprint.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_fingerprint.py flags.
    """
    p = argparse.ArgumentParser(description="Plot fingerprints.")
    p.add_argument("--chip", nargs="+", required=True, help="ChIP fingerprint TSVs.")
    p.add_argument("--control", nargs="*", default=[], help="Control fingerprint TSVs.")
    p.add_argument(
        "-o", "--output", required=True, help="SVG to write (PNG beside it)."
    )
    return p


def main(argv=None):
    """Draw one curve per library and the uniform diagonal.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    pk.use("plain")
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.plot(
        [0, 1],
        [0, 1],
        color=pk.PALETTES["plain"].colors["gray"],
        ls=":",
        label="uniform",
    )
    curves = [(p, "-", "") for p in args.chip] + [
        (p, "--", "control ") for p in args.control
    ]
    for path, style, prefix in curves:
        x, y = _curve(path)
        ax.plot(x, y, ls=style, label=prefix + Path(path).name.split(".")[1])
    ax.set(xlabel="Fraction of bins, by read count", ylabel="Fraction of reads")
    fig.legend(*ax.get_legend_handles_labels(), loc="outside right upper")
    pk.fit_legend(fig)
    pk.despine(ax)
    pk.save_figure(fig, args.output)
    plt.close(fig)


def _curve(path):
    """Turn a bin count distribution into a fingerprint curve.

    Parameters
    ----------
    path : str
        TSV with `reads` (count per bin) and `bins` (bins with that count).

    Returns
    -------
    tuple of numpy.ndarray
        Cumulative bin fraction and cumulative read fraction, from (0, 0).
        Bins with equal counts add reads at the same rate, so straight lines
        between the points are exact.
    """
    t = pl.read_csv(path, separator="\t").sort("reads")
    bins = t["bins"].to_numpy()
    reads = t["reads"].to_numpy() * bins
    x = np.concatenate([[0], np.cumsum(bins)]) / bins.sum()
    y = np.concatenate([[0], np.cumsum(reads)]) / max(reads.sum(), 1)
    return x, y


if __name__ == "__main__":
    main()
