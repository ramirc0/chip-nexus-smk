"""Gather every raw QC file of a run into one long table.

Each input is named `<exp>.<name>.<suffix>`, where `<name>` is a replicate,
`pooled`, or an IDR label. The table has one row per experiment, name, and
metric. Dedup metrics come from the `.dedup.tsv` histogram only, so any dedup
method that writes it is supported.
"""

import argparse
import json
import re
from pathlib import Path

import polars as pl
from zlib_ng import gzip_ng


def build_parser():
    """Build the qc_table.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the qc_table.py flags.
    """
    p = argparse.ArgumentParser(description="Gather QC files into one table.")
    kinds = {
        "barcode": "cutadapt pass-1 (barcode) JSON reports.",
        "adapter": "cutadapt pass-2 (adapter) JSON reports.",
        "bowtie2": "bowtie2 summaries.",
        "idxstats": "samtools idxstats of the filtered R1 reads.",
        "dedup": "Dedup group-size histograms.",
        "insert": "samtools stats of paired-end libraries.",
        "frip": "FRiP tables.",
        "peaks": "Blacklist-filtered replicate and pooled peaks.",
        "idr-status": "IDR status files.",
        "idr-peaks": "Thresholded IDR peaks.",
        "reproducibility": "Reproducibility tables.",
        "footprint": "Strand profiles around optimal summits.",
    }
    for kind, text in kinds.items():
        p.add_argument(f"--{kind}", nargs="*", default=[], help=text)
    p.add_argument("-o", "--output", required=True, help="QC TSV to write.")
    return p


def main(argv=None):
    """Parse each QC file and write the long table.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    parsers = [
        (args.barcode, ".cutadapt.barcode.json", _barcode),
        (args.adapter, ".cutadapt.adapter.json", _adapter),
        (args.bowtie2, ".bowtie2.txt", _bowtie2),
        (args.idxstats, ".filt.idxstats", _idxstats),
        (args.dedup, ".dedup.tsv", _dedup),
        (args.insert, ".samtools_stats.txt", _insert),
        (args.frip, ".frip.tsv", _frip),
        (args.peaks, ".bfilt.narrowPeak.gz", lambda p: {"peaks": _count(p)}),
        (
            args.idr_status,
            ".idr.status",
            lambda p: {"idr_status": Path(p).read_text().strip()},
        ),
        (args.idr_peaks, ".idr.narrowPeak.gz", lambda p: {"idr_peaks": _count(p)}),
    ]
    rows = []
    for paths, suffix, parse in parsers:
        for path in paths:
            exp, name = Path(path).name.removesuffix(suffix).split(".", 1)
            rows += [(exp, name, k, v) for k, v in parse(path).items()]
    for path in args.reproducibility:
        exp = Path(path).name.removesuffix(".reproducibility.tsv")
        table = pl.read_csv(path, separator="\t", infer_schema_length=0)
        rows += [(exp, "", k, v) for k, v in table.row(0, named=True).items()]

    for path in args.footprint:
        exp = Path(path).name.removesuffix(".footprint.tsv")
        rows += [
            (exp, n, k, v)
            for n, metrics in _footprint(path).items()
            for k, v in metrics.items()
        ]

    pl.DataFrame(
        [(e, n, k, None if v is None else str(v)) for e, n, k, v in rows],
        schema=["experiment", "replicate", "metric", "value"],
        orient="row",
    ).write_csv(args.output, separator="\t")


def _barcode(path):
    """Parse a cutadapt pass-1 report.

    Parameters
    ----------
    path : str
        The cutadapt JSON report.

    Returns
    -------
    dict
        Input reads, reads with a barcode, and reads whose fixed barcode had
        one substitution. A library without a barcode passes every read.
    """
    report = json.loads(Path(path).read_text())
    rescued = sum(
        sum(t["counts"][1:])
        for a in report["adapters_read1"] or []
        for t in a["five_prime_end"]["trimmed_lengths"]
    )
    return {
        "reads_raw": report["read_counts"]["input"],
        "reads_barcode": report["read_counts"]["output"],
        "reads_barcode_rescued": rescued,
    }


def _adapter(path):
    """Parse a cutadapt pass-2 report.

    Parameters
    ----------
    path : str
        The cutadapt JSON report.

    Returns
    -------
    dict
        Reads with an adapter, reads too short after trimming, and reads kept.
    """
    counts = json.loads(Path(path).read_text())["read_counts"]
    return {
        "reads_with_adapter": counts["read1_with_adapter"],
        "reads_too_short": counts["filtered"]["too_short"] or 0,
        "reads_trimmed": counts["output"],
    }


def _bowtie2(path):
    """Parse a bowtie2 summary.

    Parameters
    ----------
    path : str
        The bowtie2 stderr. Counts are pairs for paired-end input.

    Returns
    -------
    dict
        Reads given to bowtie2 and the overall alignment rate in percent.
    """
    text = Path(path).read_text()
    return {
        "reads_aligner_input": int(re.search(r"(\d+) reads; of these:", text)[1]),
        "alignment_rate": float(
            re.search(r"([\d.]+)% overall alignment rate", text)[1]
        ),
    }


def _idxstats(path):
    """Sum the filtered R1 reads over every chromosome.

    Parameters
    ----------
    path : str
        The samtools idxstats output.

    Returns
    -------
    dict
        Filtered reads, before the allow-list.
    """
    table = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "length", "mapped", "unmapped"],
    )
    return {"reads_filtered": table["mapped"].sum()}


def _dedup(path):
    """Compute the ENCODE library complexity metrics from a dedup histogram.

    Parameters
    ----------
    path : str
        TSV with `group_size` and `n_groups`.

    Returns
    -------
    dict
        Total and distinct reads, NRF, PBC1, PBC2, and the duplicate
        fraction. A ratio with a zero denominator is None.
    """
    hist = pl.read_csv(
        path,
        separator="\t",
        schema_overrides={"group_size": pl.Int64, "n_groups": pl.Int64},
    )
    total = (hist["group_size"] * hist["n_groups"]).sum()
    distinct = hist["n_groups"].sum()
    m1 = hist.filter(pl.col("group_size") == 1)["n_groups"].sum()
    m2 = hist.filter(pl.col("group_size") == 2)["n_groups"].sum()
    return {
        "reads_dedup_input": total,
        "reads_distinct": distinct,
        "nrf": distinct / total if total else None,
        "pbc1": m1 / distinct if distinct else None,
        "pbc2": m1 / m2 if m2 else None,
        "duplicate_fraction": 1 - distinct / total if total else None,
    }


def _insert(path):
    """Parse the insert size summary of samtools stats.

    Parameters
    ----------
    path : str
        The samtools stats output.

    Returns
    -------
    dict
        Mean and standard deviation of the insert size.
    """
    sn = dict(
        line.split("\t")[1:3]
        for line in Path(path).read_text().splitlines()
        if line.startswith("SN\t")
    )
    return {
        "insert_size_mean": float(sn["insert size average:"]),
        "insert_size_sd": float(sn["insert size standard deviation:"]),
    }


def _frip(path):
    """Compute FRiP for each peak set of a FRiP table.

    Parameters
    ----------
    path : str
        TSV with `peak_set`, `n_reads`, and `n_in_peaks`.

    Returns
    -------
    dict
        `frip_<peak_set>` for each row. None without reads.
    """
    table = pl.read_csv(path, separator="\t")
    return {
        f"frip_{r['peak_set']}": r["n_in_peaks"] / r["n_reads"]
        if r["n_reads"]
        else None
        for r in table.iter_rows(named=True)
    }


def _footprint(path):
    """Summarize each sample's strand profile.

    Parameters
    ----------
    path : str
        TSV with `sample`, `offset`, `plus`, and `minus`.

    Returns
    -------
    dict of str to dict
        Per sample: the offset of the + and - strand maxima, and the strand
        ratio. The ratio is reads on the expected side (+ left of the summit,
        - right of it) over reads on the other side. None without reads.
    """
    table = pl.read_csv(path, separator="\t")
    out = {}
    for sample, t in table.group_by("sample", maintain_order=True):
        left, right = t.filter(pl.col("offset") < 0), t.filter(pl.col("offset") > 0)
        expected = left["plus"].sum() + right["minus"].sum()
        other = left["minus"].sum() + right["plus"].sum()
        has_reads = expected + other > 0
        out[sample[0]] = {
            "footprint_plus_mode": t["offset"][t["plus"].arg_max()]
            if has_reads
            else None,
            "footprint_minus_mode": t["offset"][t["minus"].arg_max()]
            if has_reads
            else None,
            "footprint_strand_ratio": expected / other if other else None,
        }
    return out


def _count(path):
    """Count the lines of a gzipped file.

    Parameters
    ----------
    path : str
        Gzipped text file.

    Returns
    -------
    int
        Number of lines.
    """
    with gzip_ng.open(path, "rt") as f:
        return sum(1 for _ in f)


if __name__ == "__main__":
    main()
