"""Grade a run's QC table against the ENCODE TF ChIP-seq audits.

Cutoffs, operators, levels, and audit names copy the ENCODE portal audit code,
ENCODE-DCC/encoded at 9bd8ff1 (`src/encoded/audit/analysis.py` and
`experiment.py`). Each check writes one row: level `OK`, or the level and name
of the audit it raises. Counts depend on genome size, so they are graded only
on the assemblies the portal grades. Ratios and the design are graded on every
genome.

An experiment's `status` row holds the worst level of its rows and its
libraries' rows, except the design rows `replicates` and `control`.
"""

import argparse
import math

import polars as pl

LEVELS = ["OK", "WARNING", "NOT_COMPLIANT", "ERROR"]
ENCODE_ASSEMBLIES = {"GRCh38", "hg19", "mm10", "mm9"}
MODERN_ASSEMBLIES = {"ce10", "ce11", "dm3", "dm6"}
DESIGN = {"replicates", "control"}

# (bound, level, audit) in ascending order. The first with value < bound applies.
DEPTH = [
    (5_000_000, "ERROR", "extremely low read depth"),
    (10_000_000, "NOT_COMPLIANT", "insufficient read depth"),
    (20_000_000, "WARNING", "low read depth"),
]
CONTROL_DEPTH = [(bound, level, f"control {audit}") for bound, level, audit in DEPTH]
MODERN_DEPTH = [(3_000_000, "NOT_COMPLIANT", "insufficient read depth")]
READ_LENGTH = [
    (26, "ERROR", "extremely low read length"),
    (36, "NOT_COMPLIANT", "insufficient read length"),
    (50, "WARNING", "low read length"),
]
COMPLEXITY = {
    "nrf": [
        (0.5, "NOT_COMPLIANT", "poor library complexity"),
        (0.8, "WARNING", "moderate library complexity"),
    ],
    "pbc1": [
        (0.5, "NOT_COMPLIANT", "severe bottlenecking"),
        (0.9, "WARNING", "mild to moderate bottlenecking"),
    ],
    "pbc2": [
        (1, "NOT_COMPLIANT", "severe bottlenecking"),
        (10, "WARNING", "mild to moderate bottlenecking"),
    ],
}
REPLICATES = [(2, "NOT_COMPLIANT", "unreplicated experiment")]
# peak_sets.py labels the IDR ratios with the portal's `> 2` rule.
REPRODUCIBILITY = {
    "pass": ("OK", None),
    "borderline": ("WARNING", "borderline replicate concordance"),
    "fail": ("NOT_COMPLIANT", "insufficient replicate concordance"),
}


def build_parser():
    """Build the encode_audit.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the encode_audit.py flags.
    """
    p = argparse.ArgumentParser(description="Grade a QC table against ENCODE.")
    p.add_argument("--qc", required=True, help="Run QC TSV from qc_table.py.")
    p.add_argument("--experiments", nargs="+", required=True, help="Every experiment.")
    p.add_argument(
        "--genomes",
        nargs="+",
        required=True,
        help="Genome key of each experiment, named after its assembly.",
    )
    p.add_argument(
        "--controls",
        nargs="*",
        default=[],
        help="ChIP experiment and control pairs, flattened.",
    )
    p.add_argument("-o", "--output", required=True, help="Audit TSV to write.")
    return p


def main(argv=None):
    """Grade every library and experiment, and write the audit table.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    table = pl.read_csv(args.qc, separator="\t", infer_schema_length=0)
    qc = {
        (exp, rep or "", metric): value for exp, rep, metric, value in table.iter_rows()
    }
    reps = {}
    for exp, rep, metric in qc:
        if metric == "reads_distinct":
            reps.setdefault(exp, []).append(rep)

    controls = dict(zip(args.controls[::2], args.controls[1::2], strict=True))

    rows = []
    for exp, genome in zip(args.experiments, args.genomes, strict=True):
        graded = [
            row
            for rep in reps[exp]
            for row in _library_rows(qc, exp, rep, _library_tiers(genome))
        ]
        if exp not in controls.values():
            graded += _experiment_rows(qc, exp, genome, controls.get(exp, ""), reps)
        status = max(
            (level for _, _, metric, _, level, _ in graded if metric not in DESIGN),
            key=LEVELS.index,
            default="OK",
        )
        rows += [*graded, (exp, "", "status", None, status, None)]

    pl.DataFrame(
        rows,
        schema=dict.fromkeys(
            ["experiment", "replicate", "metric", "value", "level", "audit"],
            pl.String,
        ),
        orient="row",
    ).write_csv(args.output, separator="\t")


def _library_tiers(genome):
    """Return the tiers that grade a library.

    Parameters
    ----------
    genome : str
        Genome key, named after its assembly.

    Returns
    -------
    dict of str to list
        Each qc.tsv metric and its tiers, in output order.
    """
    if genome in ENCODE_ASSEMBLIES:
        counts = {"reads_distinct": DEPTH, "read_length_mean": READ_LENGTH}
    elif genome in MODERN_ASSEMBLIES:
        counts = {"reads_distinct": MODERN_DEPTH}
    else:
        counts = {}
    return counts | COMPLEXITY


def _library_rows(qc, exp, rep, tiers):
    """Grade one library.

    Parameters
    ----------
    qc : dict
        `(experiment, replicate, metric)` to the qc.tsv value.
    exp : str
        Experiment.
    rep : str
        Replicate.
    tiers : dict of str to list
        Metric to its tiers, from `_library_tiers`.

    Returns
    -------
    list of tuple
        One output row per metric with a value.
    """
    values = {m: qc[exp, rep, m] for m in tiers}
    numbers = {m: None if v is None else float(v) for m, v in values.items()}
    # qc_table.py writes no PBC2 when M2 = 0. With M1 > 0 it is the portal's
    # 'Infinity'.
    if numbers["pbc2"] is None and (numbers["pbc1"] or 0) > 0:
        numbers["pbc2"] = math.inf
    return [
        (exp, rep, m, values[m], *_grade(n, tiers[m]))
        for m, n in numbers.items()
        if n is not None
    ]


def _experiment_rows(qc, exp, genome, control, reps):
    """Grade the design and replicate concordance of a ChIP experiment.

    Parameters
    ----------
    qc : dict
        `(experiment, replicate, metric)` to the qc.tsv value.
    exp : str
        ChIP experiment.
    genome : str
        Its genome key.
    control : str
        Its control experiment, or empty.
    reps : dict of str to list of str
        Experiment to its replicates.

    Returns
    -------
    list of tuple
        Output rows for `replicates` and `control`. Then the control's summed
        depth on the ENCODE assemblies, and `reproducibility` with two or
        more replicates.
    """
    n = len(reps[exp])
    found = ("OK", None) if control else ("ERROR", "missing possible_controls")
    checks = [
        ("replicates", str(n), _grade(n, REPLICATES)),
        ("control", control or None, found),
    ]
    if control and genome in ENCODE_ASSEMBLIES:
        depth = sum(int(qc[control, r, "reads_distinct"]) for r in reps[control])
        checks.append(
            ("control_reads_distinct", str(depth), _grade(depth, CONTROL_DEPTH))
        )
    if n > 1:
        label = qc[exp, "", "reproducibility"]
        checks.append(("reproducibility", label, REPRODUCIBILITY[label]))
    return [(exp, "", metric, value, *grade) for metric, value, grade in checks]


def _grade(value, tiers):
    """Return the level and audit a value earns.

    Parameters
    ----------
    value : float
        Metric value.
    tiers : list of tuple
        `(bound, level, audit)` in ascending order of bound.

    Returns
    -------
    tuple of str
        Level and audit of the first tier with `value < bound`, else
        `("OK", None)`.
    """
    return next(
        ((level, audit) for bound, level, audit in tiers if value < bound),
        ("OK", None),
    )


if __name__ == "__main__":
    main()
