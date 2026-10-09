"""Write a synthetic ChIP-nexus run: genome, reads, sheet, config, and expected counts.

The genome is random, so every read aligns uniquely. Reads pile up around
planted binding sites, with + strand 5' ends left of each site and - strand
5' ends right of it. The reads cover the cases the workflow must handle:
barcodes that are exact, one substitution off, N, or junk; PCR duplicates;
read-through into the adapter; an off-list chromosome; paired-end reads; a
barcode set; libraries without a barcode; and a control.

`expected.json` holds the counts the workflow should report, derived here
from the read keys, independently of the workflow. One step is not modeled:
the random genome can resemble the adapter, so cutadapt itself decides which
reads adapter trimming leaves too short.
"""

import argparse
import json
import random
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

from zlib_ng import gzip_ng

ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
# The workflow's pass-2 adapters (preprocess.smk, common.smk).
TRIM_R1 = [ADAPTER_R1[i:22] for i in range(5)]
TRIM_R2_TAIL = ADAPTER_R2[:33]
CHROMS = {"chr1": 400_000, "chr2": 200_000, "chrOff": 50_000}
KEPT = ("chr1", "chr2")
BLACKLIST = ("chr1", 10_000, 30_000)
READ_LEN = 50
COMP = str.maketrans("ACGTN", "TGCAN")

# name -> replicates, layout, barcode patterns ("none" = no barcode), control
EXPERIMENTS = {
    "chipA": {
        "reps": {"rep1": 2, "rep2": 1},
        "layout": "se",
        "barcodes": "",
        "control": "ctlA",
    },
    "ctlA": {"reps": {"rep1": 1}, "layout": "se", "barcodes": "none", "control": ""},
    "chipB": {
        "reps": {"rep1": 1, "rep2": 1},
        "layout": "pe",
        "barcodes": "NNNNNCTGA,NNNNNTGAC,NNNNNGACT,NNNNNACTG",
        "control": "",
    },
    "chipC": {"reps": {"rep1": 1}, "layout": "se", "barcodes": "none", "control": ""},
    # Too few sites for IDR, which needs 20 merged peaks.
    "chipWeak": {
        "reps": {"rep1": 1, "rep2": 1},
        "layout": "se",
        "barcodes": "",
        "control": "",
        "sites": 8,
    },
}


def build_parser():
    """Build the make_fixture.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the make_fixture.py flags.
    """
    p = argparse.ArgumentParser(description="Write a synthetic ChIP-nexus run.")
    p.add_argument(
        "-o", "--outdir", default="resources/fixture", help="Output directory."
    )
    p.add_argument("--sites", type=int, default=200, help="Planted binding sites.")
    p.add_argument("--seed", type=int, default=0, help="Random seed.")
    return p


def main(argv=None):
    """Write the genome, reads, sheet, config, and expected counts.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    rng = random.Random(args.seed)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    genome = {c: "".join(rng.choices("ACGT", k=n)) for c, n in CHROMS.items()}
    _write_genome(genome, out)
    sites = _sites(rng, args.sites)

    controls = {e["control"] for e in EXPERIMENTS.values()}
    rows, expected = [], {}
    for exp, e in EXPERIMENTS.items():
        for rep, n_runs in e["reps"].items():
            reads = _library(
                rng, genome, [] if exp in controls else sites[: e.get("sites")], e
            )
            _mark_trimmed_short(reads, e)
            paths = _write_runs(reads, out, f"{exp}_{rep}", n_runs, e["layout"])
            for r1, r2 in paths:
                rows.append([exp, rep, r1, r2, "fixture", e["barcodes"], e["control"]])
            mode = (
                "umi"
                if e["barcodes"] != "none"
                else ("position" if exp in controls else "none")
            )
            expected[f"{exp}.{rep}"] = _expect(reads, e["layout"], mode)

    header = [
        "experiment",
        "replicate",
        "fastq_1",
        "fastq_2",
        "genome",
        "barcodes",
        "control",
    ]
    (out / "samples.tsv").write_text(
        "\n".join("\t".join(r) for r in [header, *rows]) + "\n"
    )
    (out / "config.yaml").write_text(
        f"samples: {out}/samples.tsv\n"
        "outdir: results\n"
        "pixi_env: default\n"
        "run_id: fixture\n"
        "barcodes: NNNNNCTGA\n"
        "smooth_window: 150\n"
        "genomes:\n"
        "  fixture:\n"
        f"    fasta: {out}/genome.fa\n"
        f"    chrom_sizes: {out}/genome.minimal.chrom.sizes\n"
        f"    blacklist: {out}/blacklist.bed\n"
        f"    gsize: {sum(CHROMS[c] for c in KEPT)}\n"
    )
    (out / "expected.json").write_text(json.dumps(expected, indent=1) + "\n")


def _write_genome(genome, out):
    """Write the FASTA, its index, the allow-list, and the blacklist.

    Parameters
    ----------
    genome : dict of str to str
        Sequence per chromosome.
    out : pathlib.Path
        Output directory.
    """
    fasta, fai, offset = [], [], 0
    for chrom, seq in genome.items():
        offset += len(f">{chrom}\n")
        lines = [seq[i : i + 60] for i in range(0, len(seq), 60)]
        fai.append(f"{chrom}\t{len(seq)}\t{offset}\t60\t61")
        fasta += [f">{chrom}", *lines]
        offset += sum(len(line) + 1 for line in lines)
    (out / "genome.fa").write_text("\n".join(fasta) + "\n")
    (out / "genome.fa.fai").write_text("\n".join(fai) + "\n")
    (out / "genome.minimal.chrom.sizes").write_text(
        "".join(f"{c}\t{len(genome[c])}\n" for c in KEPT)
    )
    (out / "blacklist.bed").write_text("{}\t{}\t{}\n".format(*BLACKLIST))


def _sites(rng, n):
    """Place binding sites on the kept chromosomes, 2.5 kb apart.

    Parameters
    ----------
    rng : random.Random
        Random source.
    n : int
        Number of sites.

    Returns
    -------
    list of tuple
        `(chrom, position, strength)` per site.
    """
    slots = [(c, p) for c in KEPT for p in range(5_000, CHROMS[c] - 5_000, 2_500)]
    return [(c, p, rng.randint(15, 40)) for c, p in rng.sample(slots, n)]


def _library(rng, genome, sites, e):
    """Draw the reads of one library.

    Parameters
    ----------
    rng : random.Random
        Random source.
    genome : dict of str to str
        Sequence per chromosome.
    sites : list of tuple
        Binding sites. Empty for a control.
    e : dict
        Experiment entry of `EXPERIMENTS`.

    Returns
    -------
    list of dict
        One read (or pair) per item, with its sequences, fate, and dedup key.
    """
    starts = []
    for chrom, pos, strength in sites:
        for _ in range(strength):
            strand = rng.choice("+-")
            offset = rng.randint(3, 15)
            starts.append(
                (chrom, pos - offset if strand == "+" else pos + offset, strand)
            )
    for _ in range(4000):
        chrom = rng.choices(list(CHROMS), weights=list(CHROMS.values()))[0]
        starts.append(
            (chrom, rng.randint(1_000, CHROMS[chrom] - 1_000), rng.choice("+-"))
        )

    patterns = (
        [] if e["barcodes"] == "none" else (e["barcodes"] or "NNNNNCTGA").split(",")
    )
    reads = []
    for chrom, five, strand in starts:
        read = _read(rng, genome, chrom, five, strand, patterns, e["layout"])
        reads.append(read)
        # PCR duplicates: same molecule, same UMI.
        if rng.random() < 0.15:
            reads += [dict(read) for _ in range(rng.randint(1, 3))]
        # SE duplicate with a shorter insert: still the same 5' end.
        if e["layout"] == "se" and read["fate"] == "kept" and rng.random() < 0.05:
            reads.append(
                _read(rng, genome, chrom, five, strand, patterns, "se", read, frag=30)
            )
    return reads


def _read(rng, genome, chrom, five, strand, patterns, layout, like=None, frag=None):
    """Draw one read or pair whose R1 5' end is at `five`.

    Parameters
    ----------
    rng : random.Random
        Random source.
    genome : dict of str to str
        Sequence per chromosome.
    chrom : str
        Chromosome.
    five : int
        0-based R1 5' end.
    strand : str
        `+` or `-`.
    patterns : list of str
        Barcode patterns. Empty for none.
    layout : str
        `se` or `pe`.
    like : dict, optional
        Read whose barcode and UMI to copy.
    frag : int, optional
        Fragment length. Drawn when absent.

    Returns
    -------
    dict
        `r1`, `r2` (None for SE), `fate` (`kept`, `junk`, `short`), and `key`.
    """
    fate, prefix, key_bc = "kept", "", ""
    if patterns:
        if like:
            prefix, fate, key_bc = like["prefix"], like["fate"], like["key_bc"]
        else:
            pattern = rng.choice(patterns)
            fixed = pattern.lstrip("N")
            umi = "".join(rng.choices("ACGT", k=pattern.count("N")))
            if rng.random() < 0.03:
                umi = umi[:2] + "N" + umi[3:]
            kind = rng.choices(["exact", "sub", "n", "junk"], [85, 5, 2, 8])[0]
            seen = list(fixed)
            if kind in ("sub", "junk"):
                for i in rng.sample(range(len(fixed)), 1 if kind == "sub" else 2):
                    seen[i] = rng.choice([b for b in "ACGT" if b != fixed[i]])
            if kind == "n":
                seen[rng.randrange(len(fixed))] = "N"
            prefix = umi + "".join(seen)
            fate = "junk" if kind == "junk" else "kept"
            key_bc = fixed + umi
    if frag is None and layout == "pe":
        frag = rng.randint(100, 250)
    elif frag is None:
        # 1-3 nt of adapter is below -O 4, stays untrimmed, and lowers MAPQ.
        room = READ_LEN - len(prefix)
        frag = rng.choice([n for n in [15, *range(25, 61)] if not 1 <= room - n <= 3])
    start = five if strand == "+" else five - frag + 1
    fragment = genome[chrom][start : start + frag]
    if strand == "-":
        fragment = fragment[::-1].translate(COMP)
    if fate == "kept" and frag < 22:
        fate = "short"

    r1 = (prefix + fragment + ADAPTER_R1)[:READ_LEN]
    r2 = None
    if layout == "pe":
        r2 = (
            fragment[::-1].translate(COMP) + prefix[::-1].translate(COMP) + ADAPTER_R2
        )[:READ_LEN]
    end2 = start if strand == "-" else start + frag - 1
    key = (chrom, five, strand, end2 if layout == "pe" else None, key_bc)
    return {
        "r1": r1,
        "r2": r2,
        "fate": fate,
        "key": key,
        "prefix": prefix,
        "key_bc": key_bc,
    }


def _mark_trimmed_short(reads, e):
    """Mark the reads that adapter trimming leaves under 22 nt.

    Runs cutadapt with the workflow's pass-2 settings on each read after its
    barcode, and marks the reads it drops.

    Parameters
    ----------
    reads : list of dict
        Reads from `_library`. Their `fate` is updated in place.
    e : dict
        Experiment entry of `EXPERIMENTS`.
    """
    patterns = (
        [] if e["barcodes"] == "none" else (e["barcodes"] or "NNNNNCTGA").split(",")
    )
    cut = len(patterns[0]) if patterns else 0
    flags = [f for a in TRIM_R1 for f in ("-a", a)]
    if e["layout"] == "pe":
        flags.append("--interleaved")
        for p in patterns or [""]:
            flags += ["-A", p[::-1].translate(COMP) + TRIM_R2_TAIL]
    with tempfile.TemporaryDirectory() as tmp:
        fastq, short = Path(tmp, "in.fq"), Path(tmp, "short.fq")
        with open(fastq, "w") as f:
            for i, read in enumerate(reads):
                if read["fate"] == "kept":
                    mates = [read["r1"][cut:]] + (
                        [read["r2"]] if e["layout"] == "pe" else []
                    )
                    f.writelines(f"@{i}\n{m}\n+\n{'I' * len(m)}\n" for m in mates)
        subprocess.run(
            ["cutadapt", *flags, "-m", "22", "-O", "4", "-e", "0.2"]
            + ["--too-short-output", short, "-o", "/dev/null", fastq],
            check=True,
            capture_output=True,
        )
        for line in short.read_text().splitlines()[::4]:
            reads[int(line[1:])]["fate"] = "short"


def _write_runs(reads, out, name, n_runs, layout):
    """Write a library's reads, split into runs.

    Parameters
    ----------
    reads : list of dict
        Reads from `_library`.
    out : pathlib.Path
        Output directory.
    name : str
        Library name used in file names and read IDs.
    n_runs : int
        Number of runs (FASTQ files or pairs).
    layout : str
        `se` or `pe`.

    Returns
    -------
    list of tuple of str
        `(fastq_1, fastq_2)` per run. `fastq_2` is empty for SE.
    """
    paths = []
    for run in range(n_runs):
        mates = ("1", "2") if layout == "pe" else ("1",)
        files = {m: out / f"{name}_run{run}_{m}.fastq.gz" for m in mates}
        for m, path in files.items():
            with gzip_ng.open(path, "wt") as h:
                for i, read in enumerate(reads[run::n_runs]):
                    seq = read["r1"] if m == "1" else read["r2"]
                    h.write(f"@{name}.{run}.{i} {m}:N:0\n{seq}\n+\n{'I' * len(seq)}\n")
        paths.append((str(files["1"]), str(files.get("2", ""))))
    return paths


def _expect(reads, layout, mode):
    """Derive the counts the workflow should report for one library.

    Parameters
    ----------
    reads : list of dict
        Reads from `_library`.
    layout : str
        `se` or `pe`.
    mode : str
        Dedup mode: `umi`, `position`, or `none`.

    Returns
    -------
    dict
        Barcode pass count, deduplicated reads, and track sums per strand,
        on kept chromosomes.
    """
    passed = [r for r in reads if r["fate"] != "junk"]
    kept = [r for r in passed if r["fate"] == "kept"]
    on_list = [r["key"] for r in kept if r["key"][0] in KEPT]
    if mode == "position":
        on_list = [k[:4] for k in on_list]
    groups = Counter(on_list)
    keys = on_list if mode == "none" else list(groups)
    return {
        "reads_raw": len(reads),
        "reads_barcode": len(passed),
        "reads_trimmed": len(kept),
        "reads_dedup_input": len(on_list),
        "reads_distinct": len(groups),
        "track_plus": sum(k[2] == "+" for k in keys),
        "track_minus": sum(k[2] == "-" for k in keys),
    }


if __name__ == "__main__":
    main()
