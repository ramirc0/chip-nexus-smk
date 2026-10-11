"""Check a finished run for internal consistency, and against expected counts.

For each library, the run MUST satisfy:

- cutadapt pass-2 output count == bowtie2 input count.
- Dedup histogram total == filtered R1 reads on kept chromosomes. This is the
  dedup contract, so any dedup method MUST pass it. Single-end only.
- bigWig sum == tagAlign lines == deliverable R1 reads on kept chromosomes.
  Pools sum their replicates.
- Under 0.5% of deliverable reads contain the adapter.

With `--expected` (from make_fixture.py), the QC table and track sums MUST
also match exactly. Exits non-zero if any check fails.
"""

import argparse
import json
import sys
from pathlib import Path

import polars as pl
import pyBigWig
import pysam
from zlib_ng import gzip_ng

ADAPTER = "AGATCGGAAGAGC"
ADAPTER_RC = "GCTCTTCCGATCT"


def build_parser():
    """Build the check_run.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the check_run.py flags.
    """
    p = argparse.ArgumentParser(description="Check a finished run.")
    p.add_argument("rundir", help="Run directory, results/<run_id>.")
    p.add_argument("--expected", help="expected.json from make_fixture.py.")
    return p


def main(argv=None):
    """Run every check and report failures.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    run = Path(args.rundir)
    qc = pl.read_csv(run / "qc.tsv", separator="\t", infer_schema_length=0)
    metric = {
        (r["experiment"], r["replicate"], r["metric"]): r["value"]
        for r in qc.iter_rows(named=True)
    }
    expected = json.loads(Path(args.expected).read_text()) if args.expected else {}

    failures = []

    def check(name, got, want):
        ok = got == want
        print(f"{'ok  ' if ok else 'FAIL'} {name}: {got} vs {want}")
        if not ok:
            failures.append(name)

    r1_counts = {}
    for bam in sorted(run.glob("*/*/*.bam")):
        exp, rep = bam.name.removesuffix(".bam").split(".")
        prefix = bam.with_suffix("")
        kept = _kept_chroms(prefix)
        r1, adapter, total = _scan_bam(bam, kept)
        r1_counts[exp, rep] = r1
        hist = pl.read_csv(f"{prefix}.dedup.tsv", separator="\t")
        stats = pl.read_csv(
            f"{prefix}.filt.idxstats",
            separator="\t",
            has_header=False,
            new_columns=["chrom", "length", "mapped", "unmapped"],
        )
        lib = f"{exp}.{rep}"
        check(
            f"{lib} trimmed == aligner input",
            metric[exp, rep, "reads_trimmed"],
            metric[exp, rep, "reads_aligner_input"],
        )
        if not Path(f"{prefix}.samtools_stats.txt").exists():
            check(
                f"{lib} dedup total == filtered on kept",
                (hist["group_size"] * hist["n_groups"]).sum(),
                stats.filter(pl.col("chrom").is_in(kept))["mapped"].sum(),
            )
        _check_tracks(check, lib, prefix, r1)
        check(f"{lib} adapter reads < 0.5%", adapter <= 0.005 * total, True)
        for key, want in expected.get(lib, {}).items():
            if key.startswith("track_"):
                got = _bw_sum(f"{prefix}.{'+' if key == 'track_plus' else '-'}.bw")
            else:
                got = int(metric[exp, rep, key])
            check(f"{lib} {key}", got, want)

    for pooled in sorted(run.glob("*/pooled/*.pooled.tagAlign.gz")):
        exp = pooled.name.split(".")[0]
        r1 = sum(n for (e, _), n in r1_counts.items() if e == exp)
        _check_tracks(
            check, f"{exp}.pooled", Path(str(pooled).removesuffix(".tagAlign.gz")), r1
        )

    print(f"{len(failures)} failed")
    sys.exit(1 if failures else 0)


def _kept_chroms(prefix):
    """Return the allow-list chromosomes, read from a track header.

    Parameters
    ----------
    prefix : pathlib.Path
        Library prefix.

    Returns
    -------
    list of str
        Chromosomes in the bigWig header.
    """
    bw = pyBigWig.open(f"{prefix}.+.bw")
    chroms = list(bw.chroms())
    bw.close()
    return chroms


def _scan_bam(path, kept):
    """Count R1 reads on kept chromosomes, and reads that carry the adapter.

    Parameters
    ----------
    path : pathlib.Path
        Deliverable BAM.
    kept : list of str
        Allow-list chromosomes.

    Returns
    -------
    tuple of int
        R1 reads on kept chromosomes, reads with the adapter, and all reads.
    """
    r1 = adapter = total = 0
    with pysam.AlignmentFile(path) as bam:
        for read in bam.fetch(until_eof=True):
            total += 1
            r1 += read.reference_name in kept and not read.is_read2
            seq = read.query_sequence or ""
            adapter += ADAPTER in seq or ADAPTER_RC in seq
    return r1, adapter, total


def _check_tracks(check, name, prefix, r1):
    """Check that track sums, tagAlign lines, and R1 reads agree.

    Parameters
    ----------
    check : callable
        Records `check(name, got, want)`.
    name : str
        Library or pool name for messages.
    prefix : pathlib.Path
        Output prefix of the library or pool.
    r1 : int
        Deliverable R1 reads on kept chromosomes.
    """
    track = _bw_sum(f"{prefix}.+.bw") + _bw_sum(f"{prefix}.-.bw")
    with gzip_ng.open(f"{prefix}.tagAlign.gz", "rt") as f:
        lines = sum(1 for _ in f)
    check(f"{name} track sum == tagAlign lines", track, lines)
    check(f"{name} tagAlign lines == R1 on kept", lines, r1)


def _bw_sum(path):
    """Sum every value in a bigWig.

    Parameters
    ----------
    path : str
        The bigWig file.

    Returns
    -------
    int
        Sum over all chromosomes.
    """
    bw = pyBigWig.open(path)
    total = sum(v for c in bw.chroms() for _, _, v in bw.intervals(c) or [])
    bw.close()
    return int(total)


if __name__ == "__main__":
    main()
