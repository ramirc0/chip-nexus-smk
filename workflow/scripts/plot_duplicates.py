"""Plot how reads split over duplicate group sizes, one line per replicate.

A complex library keeps most reads in groups of size 1. Reads piling into
large groups mean the library was sequenced past its complexity.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import poikilos as pk
import polars as pl


def build_parser():
    """Build the plot_duplicates.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_duplicates.py flags.
    """
    p = argparse.ArgumentParser(description="Plot duplicate group sizes.")
    p.add_argument(
        "-i", "--input", nargs="+", required=True, help="Dedup histogram TSVs."
    )
    p.add_argument(
        "-o", "--output", required=True, help="SVG to write (PNG beside it)."
    )
    return p


def main(argv=None):
    """Draw the fraction of reads at each group size.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    pk.use("plain")
    fig, ax = plt.subplots(figsize=(4, 3))
    for path in args.input:
        hist = pl.read_csv(path, separator="\t").with_columns(
            reads=pl.col("group_size") * pl.col("n_groups")
        )
        ax.plot(
            hist["group_size"],
            hist["reads"] / hist["reads"].sum(),
            marker="o",
            label=Path(path).name.split(".")[1],
        )
    ax.set(
        xscale="log",
        xlabel="Duplicate group size",
        ylabel="Fraction of reads",
    )
    fig.legend(*ax.get_legend_handles_labels(), loc="outside right upper")
    pk.fit_legend(fig)
    pk.despine(ax)
    pk.save_figure(fig, args.output)
    plt.close(fig)


if __name__ == "__main__":
    main()
