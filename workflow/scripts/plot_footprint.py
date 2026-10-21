"""Plot the strand footprint around summits, one panel per sample.

+ strand 5' ends are drawn up and - strand 5' ends down, the usual ChIP-nexus
view. Good data puts the + peak left of the summit and the - peak right of it.
Pools are deeper than replicates, so each panel has its own y axis.
"""

import argparse

import matplotlib.pyplot as plt
import poikilos as pk
import polars as pl


def build_parser():
    """Build the plot_footprint.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_footprint.py flags.
    """
    p = argparse.ArgumentParser(description="Plot strand footprints.")
    p.add_argument("-i", "--input", required=True, help="Footprint TSV.")
    p.add_argument(
        "-o", "--output", required=True, help="SVG to write (PNG beside it)."
    )
    return p


def main(argv=None):
    """Draw both strands for each sample.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    table = pl.read_csv(args.input, separator="\t")
    samples = table["sample"].unique(maintain_order=True).to_list()
    colors = pk.PALETTES["plain"].colors

    pk.use("plain")
    fig, axes = plt.subplots(
        1, len(samples), figsize=(3.5 * len(samples), 3), squeeze=False
    )
    for ax, sample in zip(axes[0], samples):
        t = table.filter(pl.col("sample") == sample)
        ax.plot(t["offset"], t["plus"], color=colors["blue"], label="+ strand")
        ax.plot(t["offset"], -t["minus"], color=colors["red"], label="- strand")
        pk.zero_line(ax)
        ax.set(
            title=sample,
            xlabel="Distance from summit (bp)",
            ylabel="Mean 5' ends per summit",
        )
    fig.legend(*axes[0][0].get_legend_handles_labels(), loc="outside right upper")
    pk.fit_legend(fig)
    for ax in axes[0]:
        pk.despine(ax)
    pk.save_figure(fig, args.output)
    plt.close(fig)


if __name__ == "__main__":
    main()
