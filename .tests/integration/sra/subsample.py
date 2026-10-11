"""Write the first reads of each run in the sheet, taken from the sra-chip store.

Each FASTQ path in the sheet names its output. The source is the file of the
same name under `--source`.
"""

import argparse
import csv
from itertools import islice
from pathlib import Path

from zlib_ng import gzip_ng


def build_parser():
    """Build the subsample.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the subsample.py flags.
    """
    p = argparse.ArgumentParser(description="Subsample the sheet's SRA runs.")
    p.add_argument(
        "--sheet",
        default=".tests/integration/sra/samples.tsv",
        help="Sample sheet whose FASTQ paths to write.",
    )
    p.add_argument(
        "--source",
        default="/data/zusers/ramirezc/data/sra-chip/reads",
        help="Directory with the full FASTQs.",
    )
    p.add_argument("--reads", type=int, default=2_000_000, help="Reads per file.")
    return p


def main(argv=None):
    """Write the first `--reads` records of every FASTQ in the sheet.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.sheet) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    for path in [Path(r[c]) for r in rows for c in ("fastq_1", "fastq_2") if r[c]]:
        path.parent.mkdir(parents=True, exist_ok=True)
        with (
            gzip_ng.open(Path(args.source, path.name), "rt") as src,
            gzip_ng.open(path, "wt") as out,
        ):
            out.writelines(islice(src, 4 * args.reads))


if __name__ == "__main__":
    main()
