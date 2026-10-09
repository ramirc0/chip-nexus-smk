# chip-nexus-smk

ChIP-nexus processing in the ENCODE TF ChIP-seq style. It takes FASTQs and
writes deduplicated BAMs, strand-specific 5' count tracks, MACS3 peaks, IDR
reproducibility sets, and one QC table.

The read processing follows He, Johnston et al. 2015 (*Nat. Biotechnol.*
33:395). Replicates, pseudoreplicates, IDR, and the QC metrics follow the
ENCODE pipeline.

## Steps

1. **Barcode** (cutadapt pass 1). Cut the UMI, then require the fixed barcode
   at the new read start. The read ID gains `_<KEY>`, where KEY is the
   canonical fixed barcode plus the UMI as read.
2. **Adapter** (cutadapt pass 2). Trim the 3' adapter and keep reads of 22 nt
   or more. Paired-end R2 also loses the reverse-complemented barcode.
3. **Align** with bowtie2 end-to-end. Keep MAPQ >= 30 and `-F 1804`. Pairs must
   be proper, and orphaned mates are dropped.
4. **Dedup** with `samtools markdup`, keyed on the read ID UMI. The BAM keeps
   R1 only.
5. **Tracks** with bam2bw: raw 5' end counts, one bigWig per strand.
6. **Peaks** with MACS3 on each replicate, the pool, and their pseudoreplicates.
7. **IDR** on true replicate pairs, self-pseudoreplicates, and pooled
   pseudoreplicates. Then pick the ENCODE conservative and optimal sets.
8. **QC**: per-step read counts, NRF, PBC1, PBC2, FRiP, and reproducibility.

## Sample sheet

One row per sequencing run, tab-separated. Runs of the same experiment and
replicate are concatenated.

| Column | Required | Meaning |
|---|---|---|
| `experiment` | yes | Experiment name. |
| `replicate` | yes | Replicate name. `pooled`, `peaks`, `idr`, and `pr` are reserved. |
| `fastq_1` | yes | Read 1 FASTQ, gzipped or not. |
| `fastq_2` | no | Read 2 FASTQ. Makes the library paired-end. |
| `genome` | yes | Key in `config["genomes"]`. |
| `barcodes` | no | 5' read pattern, N = UMI base, e.g. `NNNNNCTGA`. Comma-separate a barcode set. Empty uses `config["barcodes"]`. `none` means no barcode. |
| `control` | no | Experiment to use as the MACS3 control. |

Names match `[A-Za-z0-9][A-Za-z0-9_-]*`. Dots separate name parts in output
file names, so names cannot contain one.

A control experiment gets BAMs and tracks but no peaks. MACS3 pools all of its
replicates.

## Config

`config/config.yaml.template` lists the five supported genomes (GRCh38, mm10,
dm6, sacCer3, and ASM294v2). Each needs:

- `fasta`, with a `samtools faidx` index beside it.
- `chrom_sizes`, the allow-list. Everything after dedup keeps only these
  chromosomes. The workflow checks names and lengths against the `.fai`.
- `blacklist`, applied to peaks only. Use an empty file for none.
- `gsize`, the MACS3 effective genome size.

`smooth_window` (150) sets the MACS3 window centered on each 5' end.

## Outputs

Under `results/<run_id>/<experiment>/`:

| Path | What |
|---|---|
| `<rep>/<exp>.<rep>.bam` | Filtered, deduplicated R1 reads, all chromosomes. |
| `<rep or pooled>/<exp>.<sample>.{+,-}.bw` | 5' end counts per strand, both positive. |
| `<rep>/<exp>.<rep>.cutadapt.{barcode,adapter}.json` | cutadapt reports. |
| `<rep>/<exp>.<rep>.bowtie2.txt` | bowtie2 summary. |
| `<rep>/<exp>.<rep>.dedup.tsv` | Duplicate group-size histogram. |
| `<rep>/<exp>.<rep>.markdup.json` | samtools markdup stats. |
| `peaks/<exp>.<set>.bfilt.narrowPeak.gz` | MACS3 peaks after the blacklist. |
| `idr/<exp>.<pair>.idr.narrowPeak.gz` | IDR peaks at 0.05, after the blacklist. |
| `<exp>.{conservative,optimal}.narrowPeak.gz` | ENCODE reproducible peak sets. |
| `<exp>.reproducibility.tsv` | Peak counts, rescue and self-consistency ratios. |

`results/<run_id>/qc.tsv` holds every metric in long form: `experiment`,
`replicate`, `metric`, `value`.

## Choices

| Choice | Reason |
|---|---|
| Fixed barcode allows 1 substitution, and N in the read matches anything. | Near misses are substitutions, 1.5-13% of reads. |
| R1 adapter is the 22 nt prefix `AGATCGGAAGAGCACACGTCTG`. | Reads carry the TruSeq tail, not the paper's BamHI tail. The full adapter misses read-through past about 32 nt. |
| The read's last base is kept. | The paper removes positions 1-9 only. The 2014 code's extra trim is an off-by-one. |
| Dedup key is 5' end, strand, and UMI (SE), or both ends, orientation, and UMI (PE). | The paper's key. Unlike the 2014 code, read length is not part of it. |
| A ChIP library with no barcode is not deduplicated. A control without one is deduplicated by position. | Without a UMI, position dedup caps signal at peaks. |
| MACS3 `--nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 0.01`. | Centers a 150 bp window on each 5' end. Dedup already happened. |
| Peaks capped at 500,000, IDR 0.05 ranked by p-value. | ENCODE TF ChIP-seq settings. |
| Tracks hold R1 only, as raw counts. | The R1 5' end is the λ-exonuclease stop. R2 marks the sonication end. |
| No strand cross-correlation, no MACS signal tracks. | The 5' count tracks are the signal. |

## Running

```bash
cp config/config.yaml.template config/config.yaml  # then edit
cp config/samples.example.tsv config/samples.tsv   # then edit
pixi shell                                         # launcher env
snakemake -n -p --profile profiles/local           # dry run
snakemake --profile profiles/local                 # local
snakemake --profile profiles/slurm                 # SLURM
```

The first run builds a bowtie2 index per genome under `results/refs/`. GRCh38
takes hours; build it once and reuse `results/refs/` across runs.

Tune `set-resources` in the profiles from the `benchmarks/<run_id>/` TSVs after
a first run.

## Testing

`tests/make_fixture.py` writes a synthetic run: a random genome, reads around
200 planted sites, and the counts the workflow should report.
`tests/check_run.py` checks any finished run for internal consistency, and a
fixture run against those counts.

```bash
pixi run --manifest-path workflow/envs/pixi.toml python tests/make_fixture.py
snakemake --profile profiles/local --configfile resources/fixture/config.yaml
pixi run --manifest-path workflow/envs/pixi.toml \
    python tests/check_run.py results/fixture --expected resources/fixture/expected.json
```

## Conventions

The layout and idioms come from `snakemake-template`. Follow the
[Nextstrain Snakemake style guide][style-guide] and
[Snakemake's best-practices guide][best-practices]:

- No `run:` blocks. Custom logic lives in `workflow/scripts/`, called from `shell:`.
- Raw, triple-quoted `shell` blocks (`r"""`), one option per line.
- Every interpolation quoted with `:q`. Multi-value params passed as lists.
- Config mapped into `shell` via `params:`, never interpolated directly.
- `input:`/`output:` are literal path strings. All paths are relative.
- Logs tee'd to file and terminal (`exec &> >(tee {log:q})`).
- Every rule has `log:`, `benchmark:`, and `software:`. No `message:`.
- Scripts expose `build_parser()` and `main()`. Docstrings follow numpydoc.
  `pixi run lint` checks the scripts with ruff and `numpydoc lint`. Neither
  parses `.smk` files, so keep the rule helpers in the same style by hand.

[style-guide]: https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html
[best-practices]: https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html

## Environment

The launcher env (`pixi.toml` at the root) has Snakemake, the pixi
software-deployment plugin, and the SLURM executor plugins. The rule envs are
a second pixi workspace, `workflow/envs/`, with two environments:

- `default`: cutadapt, bowtie2, samtools, bedtools, MACS3, bam2bw, gawk, and
  polars.
- `idr`: IDR 2.0.4.2, which needs python 3.9 and numpy < 1.20.

Jobs install the envs on first use. After a manifest edit, run `pixi lock` in
`workflow/envs/`. Any change to a manifest or lock reruns every rule.
