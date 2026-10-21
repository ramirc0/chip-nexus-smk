"""Plot paired-end insert sizes, one line per replicate.

In ChIP-nexus the insert runs from the exonuclease stop (R1 5' end) to the
sonication break (R2 5' end), so it reflects digestion and shearing together.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import poikilos as pk


def build_parser():
    """Build the plot_insert_size.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_insert_size.py flags.
    """
    p = argparse.ArgumentParser(description="Plot insert sizes.")
    p.add_argument(
        "-i", "--input", nargs="+", required=True, help="samtools stats files."
    )
    p.add_argument(
        "-o", "--output", required=True, help="SVG to write (PNG beside it)."
    )
    return p


def main(argv=None):
    """Draw the insert size distribution up to its 99th percentile.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    pk.use("plain")
    fig, ax = plt.subplots(figsize=(4, 3))
    for path in args.input:
        size, pairs = _insert_sizes(path)
        fraction = pairs / max(pairs.sum(), 1)
        last = min(np.searchsorted(np.cumsum(fraction), 0.99), len(size) - 1)
        keep = size <= size[last]
        ax.plot(
            size[keep],
            fraction[keep],
            label=Path(path).name.split(".")[1],
        )
    ax.set(xlabel="Insert size (bp)", ylabel="Fraction of pairs")
    fig.legend(*ax.get_legend_handles_labels(), loc="outside right upper")
    pk.fit_legend(fig)
    pk.despine(ax)
    pk.save_figure(fig, args.output)
    plt.close(fig)


def _insert_sizes(path):
    """Read the insert size histogram from samtools stats.

    Parameters
    ----------
    path : str
        The samtools stats output. Its `IS` lines hold size and pair count.

    Returns
    -------
    tuple of numpy.ndarray
        Insert sizes and the number of pairs at each.
    """
    with open(path) as f:
        rows = [line.split("\t")[1:3] for line in f if line.startswith("IS\t")]
    table = np.array(rows, dtype=np.int64).reshape(-1, 2)
    return table[:, 0], table[:, 1]


if __name__ == "__main__":
    main()
