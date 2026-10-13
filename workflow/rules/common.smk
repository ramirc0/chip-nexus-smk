# Config, sample sheet validation, and shared helpers.

import re
from itertools import combinations
from pathlib import Path

import polars as pl
from snakemake.exceptions import WorkflowError


configfile: "config/config.yaml"


RESULTS = config["outdir"]
RUN_ID = config.get("run_id", "default")
OUTDIR = f"{RESULTS}/{RUN_ID}"          # per-run
REFDIR = f"{RESULTS}/refs"              # run-independent
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"
GENOMES = config["genomes"]

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)
IDR_ENV = pixi(workspace="../envs", env="idr", locked=True)
QC_ENV = pixi(workspace="../envs", env="qc", locked=True)

# The pixi plugin installs the envs inside each job. Parallel SLURM jobs then race
# on a fresh or relocked env, so install them once before any job starts.
onstart:
    shell("pixi install --manifest-path {workflow.basedir:q}/envs/pixi.toml --all --locked")

# Dots separate name parts in file names, so names never contain one.
NAME = r"[A-Za-z0-9][A-Za-z0-9_-]*"
# Path components and IDR labels that a replicate name would collide with.
RESERVED = {"pooled", "peaks", "idr", "pr"}


_SHEET = Path(config["samples"])
if not _SHEET.exists():
    raise WorkflowError(
        f"Sample sheet not found: {_SHEET}\n"
        "Columns: experiment, replicate, fastq_1, genome, "
        "[fastq_2], [barcodes], [control]; see config/samples.example.tsv."
    )


def _fail(message):
    """Raise a sample sheet error.

    Parameters
    ----------
    message : str
        What is wrong.

    Raises
    ------
    WorkflowError
        Always, prefixed with the sheet path.
    """
    raise WorkflowError(f"{_SHEET}: {message}")


# All-string read; optional columns default to "".
_sheet = pl.read_csv(_SHEET, separator="\t", infer_schema_length=0, comment_prefix="#")
_missing = {"experiment", "replicate", "fastq_1", "genome"} - set(_sheet.columns)
if _missing:
    _fail(f"missing columns {sorted(_missing)}")
_sheet = _sheet.with_columns(
    pl.lit("").alias(c) for c in ("fastq_2", "barcodes", "control") if c not in _sheet.columns
).with_columns(pl.all().fill_null("").str.strip_chars())
if _sheet.is_empty():
    _fail("no rows")


def _patterns(cell):
    """Parse a `barcodes` cell into 5' read patterns.

    Parameters
    ----------
    cell : str
        Comma-separated patterns, `none`, or empty for `config["barcodes"]`.

    Returns
    -------
    list of str
        Patterns such as `NNNNNCTGA`. Empty for `none`.
    """
    if cell == "none":
        return []
    patterns = (cell or config["barcodes"]).split(",")
    shapes = set()
    for p in patterns:
        m = re.fullmatch(r"(N+)([ACGT]+)", p)
        if not m:
            _fail(f"barcode pattern {p!r} is not N+ (UMI) then [ACGT]+ (fixed)")
        shapes.add((len(m[1]), len(m[2])))
    if len(shapes) > 1 or len(set(patterns)) < len(patterns):
        _fail(f"barcode set {cell!r} mixes UMI or fixed lengths, or repeats a pattern")
    return patterns


EXPERIMENTS = {}  # exp -> {"genome", "control", "reps"}
LIBRARIES = {}  # (exp, rep) -> {"runs", "layout", "patterns"}
for row in _sheet.iter_rows(named=True):
    exp, rep = row["experiment"], row["replicate"]
    for name in (exp, rep):
        if not re.fullmatch(NAME, name) or name in RESERVED:
            _fail(f"name {name!r} must match {NAME} and not be one of {sorted(RESERVED)}")
    if row["genome"] not in GENOMES:
        _fail(f"{exp}: genome {row['genome']!r} is not in config['genomes']")
    if not row["fastq_1"]:
        _fail(f"{exp} {rep}: empty fastq_1")
    for fastq in (row["fastq_1"], row["fastq_2"]):
        if fastq and not fastq.endswith(".gz"):
            _fail(f"{exp} {rep}: {fastq} must be gzipped (.gz)")
    e = EXPERIMENTS.setdefault(
        exp, {"genome": row["genome"], "control": row["control"] or None, "reps": []}
    )
    if (e["genome"], e["control"]) != (row["genome"], row["control"] or None):
        _fail(f"{exp}: rows disagree on genome or control")
    layout = "pe" if row["fastq_2"] else "se"
    patterns = _patterns(row["barcodes"])
    lib = LIBRARIES.get((exp, rep))
    if lib is None:
        e["reps"].append(rep)
        lib = LIBRARIES[exp, rep] = {"runs": [], "layout": layout, "patterns": patterns}
    if (lib["layout"], lib["patterns"]) != (layout, patterns):
        _fail(f"{exp} {rep}: runs disagree on barcodes or single/paired end")
    lib["runs"].append((row["fastq_1"], row["fastq_2"]))

CONTROLS = sorted({e["control"] for e in EXPERIMENTS.values() if e["control"]})
CHIPS = [exp for exp in EXPERIMENTS if exp not in CONTROLS]
for exp, e in EXPERIMENTS.items():
    ctl = EXPERIMENTS.get(e["control"]) if e["control"] else None
    if e["control"] and (ctl is None or ctl["genome"] != e["genome"] or ctl["control"]):
        _fail(f"{exp}: control {e['control']!r} must be an experiment with the same genome and no control")


def _check_allow_list(genome):
    """Check that a genome's allow-list matches its FASTA index.

    Parameters
    ----------
    genome : str
        Key in `config["genomes"]`.

    Raises
    ------
    WorkflowError
        If a chromosome is absent from `<fasta>.fai` or has another length.
    """
    g = GENOMES[genome]
    fai = Path(f"{g['fasta']}.fai")
    if not fai.exists():
        raise WorkflowError(f"{genome}: no FASTA index {fai}; run samtools faidx")
    lengths = dict(line.split()[:2] for line in fai.read_text().splitlines())
    for line in Path(g["chrom_sizes"]).read_text().splitlines():
        chrom, length = line.split()
        if lengths.get(chrom) != length:
            raise WorkflowError(f"{genome}: {chrom} {length} in {g['chrom_sizes']} is not in {fai}")


for _genome in {e["genome"] for e in EXPERIMENTS.values()}:
    _check_allow_list(_genome)


wildcard_constraints:
    genome=NAME,
    exp=NAME,
    rep=rf"(?!pooled(?![A-Za-z0-9_-])){NAME}",
    sample=NAME,
    peakset=rf"{NAME}(\.pr[12])?",
    pair=rf"{NAME}\.{NAME}",


def lib_prefix(exp, sample):
    """Return the output prefix of a replicate or a pool.

    Parameters
    ----------
    exp : str
        Experiment.
    sample : str
        Replicate, or `pooled`.

    Returns
    -------
    str
        Path `<outdir>/<run_id>/<exp>/<sample>/<exp>.<sample>`.
    """
    return f"{OUTDIR}/{exp}/{sample}/{exp}.{sample}"


def reps_of(exp, sample):
    """Return the replicates that make up a replicate or a pool.

    Parameters
    ----------
    exp : str
        Experiment.
    sample : str
        Replicate, or `pooled`.

    Returns
    -------
    list of str
        Every replicate for `pooled`, else `[sample]`.
    """
    return EXPERIMENTS[exp]["reps"] if sample == "pooled" else [sample]


def samples_of(exp):
    """Return the replicates of an experiment, plus `pooled` when it has several.

    Parameters
    ----------
    exp : str
        Experiment.

    Returns
    -------
    list of str
        Replicates in sheet order, then `pooled` if there are two or more.
    """
    reps = EXPERIMENTS[exp]["reps"]
    return reps + ["pooled"] if len(reps) > 1 else reps


def genome_file(exp, key):
    """Return a genome file of an experiment.

    Parameters
    ----------
    exp : str
        Experiment.
    key : str
        `fasta`, `chrom_sizes`, or `blacklist`.

    Returns
    -------
    str
        Path from `config["genomes"]`.
    """
    return GENOMES[EXPERIMENTS[exp]["genome"]][key]


def bowtie2_index(exp):
    """Return the bowtie2 index prefix of an experiment's genome.

    Parameters
    ----------
    exp : str
        Experiment.

    Returns
    -------
    str
        Path `<outdir>/refs/<genome>/<genome>`.
    """
    genome = EXPERIMENTS[exp]["genome"]
    return f"{REFDIR}/{genome}/{genome}"


def dedup_mode(exp, rep):
    """Return how a library is deduplicated.

    Parameters
    ----------
    exp : str
        Experiment.
    rep : str
        Replicate.

    Returns
    -------
    str
        `umi` with a barcode, `position` for a control without one, and
        `none` for a ChIP library without one.
    """
    if LIBRARIES[exp, rep]["patterns"]:
        return "umi"
    return "position" if exp in CONTROLS else "none"


def barcode_flags(exp, rep):
    """Return the cutadapt pass-1 flags that check and remove the barcode.

    `-u` cuts the UMI before adapter matching, so the fixed barcode is
    anchored at the new read start. The rename appends `_<KEY>` to every read
    ID, where KEY is the canonical fixed barcode plus the UMI as read.

    Parameters
    ----------
    exp : str
        Experiment.
    rep : str
        Replicate.

    Returns
    -------
    list of str
        Flag tokens. Empty for a library without a barcode, so pass 1 copies
        the reads unchanged. Interpolate them with `:q`.
    """
    lib = LIBRARIES[exp, rep]
    if not lib["patterns"]:
        return []
    umi_len = lib["patterns"][0].count("N")
    flags = ["-u", str(umi_len), "--match-read-wildcards", "--discard-untrimmed"]
    for p in lib["patterns"]:
        fixed = p[umi_len:]
        flags += ["-g", f"{fixed}=^{fixed};e=1;noindels"]
    r1 = "r1." if lib["layout"] == "pe" else ""
    return flags + ["--rename", f"{{id}}_{{{r1}adapter_name}}{{{r1}cut_prefix}}"]


def _revcomp(seq):
    """Reverse-complement a DNA sequence.

    Parameters
    ----------
    seq : str
        Bases from `ACGTN`.

    Returns
    -------
    str
        Reverse complement.
    """
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def r2_adapter_flags(exp, rep):
    """Return the cutadapt pass-2 flags for read 2.

    Read 2 runs through the reverse complement of the read 1 barcode, then
    the TruSeq read 1 primer site.

    Parameters
    ----------
    exp : str
        Experiment.
    rep : str
        Replicate.

    Returns
    -------
    list of str
        `--interleaved` and one `-A` per barcode pattern for a paired-end
        library. Empty for single-end. Interpolate them with `:q`.
    """
    lib = LIBRARIES[exp, rep]
    if lib["layout"] == "se":
        return []
    tail = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    flags = ["--interleaved"]
    for p in lib["patterns"] or [""]:
        flags += ["-A", _revcomp(p) + tail]
    return flags


def control_flags(exp):
    """Return the MACS3 control flags of a ChIP experiment.

    Parameters
    ----------
    exp : str
        ChIP experiment.

    Returns
    -------
    list of str
        `-c` and every control replicate tagAlign. MACS3 pools them. Empty
        without a control.
    """
    tagaligns = control_files(exp, "tagAlign.gz")
    return ["-c", *tagaligns] if tagaligns else []


def control_files(exp, suffix):
    """Return one file per control replicate of a ChIP experiment.

    Parameters
    ----------
    exp : str
        ChIP experiment.
    suffix : str
        File suffix after the library prefix, such as `tagAlign.gz`.

    Returns
    -------
    list of str
        `<prefix>.<suffix>` per control replicate. Empty without a control.
    """
    ctl = EXPERIMENTS[exp]["control"]
    if not ctl:
        return []
    return [f"{lib_prefix(ctl, r)}.{suffix}" for r in EXPERIMENTS[ctl]["reps"]]


def pe_reps(exp):
    """Return the paired-end replicates of an experiment.

    Parameters
    ----------
    exp : str
        Experiment.

    Returns
    -------
    list of str
        Replicates whose library is paired-end, in sheet order.
    """
    return [r for r in EXPERIMENTS[exp]["reps"] if LIBRARIES[exp, r]["layout"] == "pe"]


def peaksets(exp):
    """Return every MACS3 peak set of a ChIP experiment.

    Parameters
    ----------
    exp : str
        ChIP experiment.

    Returns
    -------
    list of str
        Each sample from `samples_of`, plus its `.pr1` and `.pr2` pseudoreps.
    """
    return [f"{s}{pr}" for s in samples_of(exp) for pr in ("", ".pr1", ".pr2")]


def idr_pairs(exp):
    """Return the IDR comparisons of a ChIP experiment.

    True replicate pairs and pooled pseudoreps use the pooled peaks as the
    oracle. Self-pseudoreps use the replicate's own peaks.

    Parameters
    ----------
    exp : str
        ChIP experiment.

    Returns
    -------
    dict of str to tuple of str
        Label to `(peakset_1, peakset_2, oracle)`. Labels are `<a>.<b>` for
        true pairs, `<rep>.pr` for self-pseudoreps, and `pooled.pr`.
    """
    reps = EXPERIMENTS[exp]["reps"]
    pairs = {f"{r}.pr": (f"{r}.pr1", f"{r}.pr2", r) for r in reps}
    if len(reps) > 1:
        pairs |= {f"{a}.{b}": (a, b, "pooled") for a, b in combinations(reps, 2)}
        pairs["pooled.pr"] = ("pooled.pr1", "pooled.pr2", "pooled")
    return pairs


def idr_peaks(exp, label):
    """Return the thresholded IDR peaks of one comparison.

    Parameters
    ----------
    exp : str
        ChIP experiment.
    label : str
        Key of `idr_pairs`.

    Returns
    -------
    str
        Path of the blacklist-filtered IDR narrowPeak.
    """
    return f"{OUTDIR}/{exp}/idr/{exp}.{label}.idr.narrowPeak.gz"
