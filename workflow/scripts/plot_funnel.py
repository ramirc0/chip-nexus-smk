"""Plot where an experiment's reads go, one stacked bar per replicate.

Each bar splits the raw reads (pairs for paired-end) by the step that removed
them, ending with the distinct reads. For a library without a barcode, the
duplicate segment counts position duplicates, which a ChIP library keeps.
"""

import argparse

import matplotlib.pyplot as plt
import poikilos as pk
import polars as pl

# Segment label, then the counts before and after the step, from qc.tsv.
STEPS = [
    ("No barcode", "reads_raw", "reads_barcode"),
    ("Too short after trimming", "reads_barcode", "reads_trimmed"),
    ("Unaligned or filtered", "reads_trimmed", "reads_filtered"),
    ("Off allow-list", "reads_filtered", "reads_dedup_input"),
    ("Duplicate", "reads_dedup_input", "reads_distinct"),
    ("Distinct", "reads_distinct", None),
]


def build_parser():
    """Build the plot_funnel.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_funnel.py flags.
    """
    p = argparse.ArgumentParser(description="Plot the read funnel.")
    p.add_argument("-i", "--input", required=True, help="Run QC TSV.")
    p.add_argument("-e", "--experiment", required=True, help="Experiment to plot.")
    p.add_argument(
        "-o", "--output", required=True, help="SVG to write (PNG beside it)."
    )
    return p


def main(argv=None):
    """Draw one stacked bar per replicate of the experiment.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    counts = (
        pl.read_csv(args.input, separator="\t", infer_schema_length=0)
        .filter(
            pl.col("experiment") == args.experiment,
            pl.col("metric").is_in([s[1] for s in STEPS]),
        )
        .with_columns(pl.col("value").cast(pl.Float64))
        .pivot("metric", index=["experiment", "replicate"], values="value")
    )
    names = counts["replicate"].to_list()

    pk.use("plain")
    fig, ax = plt.subplots(figsize=(6, 1.5 + 0.3 * len(names)))
    left = 0
    for label, before, after in STEPS:
        lost = counts[before] - (counts[after] if after else 0)
        width = lost / counts["reads_raw"]
        ax.barh(names, width, left=left, label=label)
        left = left + width
    ax.invert_yaxis()
    ax.set(xlabel="Fraction of raw reads")
    fig.legend(*ax.get_legend_handles_labels(), loc="outside right upper")
    pk.fit_legend(fig)
    pk.despine(ax, categorical_y=True)
    pk.save_figure(fig, args.output)
    plt.close(fig)


if __name__ == "__main__":
    main()
