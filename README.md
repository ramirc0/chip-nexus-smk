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
5. **Tracks** with `figwig bam2bw`: raw 5' end counts, one bigWig per strand.
6. **Peaks** with MACS3 on each replicate, the pool, and their pseudoreplicates.
7. **IDR** on true replicate pairs, self-pseudoreplicates, and pooled
   pseudoreplicates. Then pick the ENCODE conservative and optimal sets.
8. **QC**: per-step read counts, read length, NRF, PBC1, PBC2, FRiP,
   reproducibility, and the strand footprint at summits. Plus the figures
   listed under Outputs.
9. **ENCODE audits**: grade each library and experiment against the ENCODE TF
   ChIP-seq audits.

## Sample sheet

One row per sequencing run, tab-separated. Runs of the same experiment and
replicate are concatenated.

| Column | Required | Meaning |
|---|---|---|
| `experiment` | yes | Experiment name. |
| `replicate` | yes | Replicate name. `pooled`, `peaks`, `idr`, and `pr` are reserved. |
| `fastq_1` | yes | Read 1 FASTQ, gzipped. |
| `fastq_2` | no | Read 2 FASTQ, gzipped. Makes the library paired-end. |
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

`run_id` names the output directory, `results/<run_id>/`. Production runs use
`<YYYY-MM-DD>-<slug>`, the start date plus a short name, such as
`2026-10-11-sra-chip`. Rerun under the same `run_id` when the sample sheet
grows. Snakemake then skips the finished libraries.

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
| `<exp>.footprint.{tsv,svg}` | Mean + and - strand 5' ends within 100 bp of the top 1,000 optimal summits. + ends should peak left of the summit and - ends right of it. |
| `<exp>.fingerprint.svg` | Cumulative reads over 500 bp bins sorted by count, ChIP and control. Data in `<rep>/<exp>.<rep>.fingerprint.tsv`. |
| `<exp>.duplicates.svg` | Fraction of reads per duplicate group size, per replicate. |
| `<exp>.funnel.svg` | Where each replicate's reads go, from raw reads to distinct reads. |
| `<exp>.insert_size.svg` | Insert sizes of paired-end replicates. |
| `idr/<exp>.<pair>.idr.svg` | The four `idr --plot` panels, redrawn from the IDR output. |

`results/<run_id>/qc.tsv` holds every metric in long form: `experiment`,
`replicate`, `metric`, `value`. `results/<run_id>/encode.tsv` grades those
metrics, as described under ENCODE audits.

Every figure is an SVG with a PNG beside it.

## ENCODE audits

`encode.tsv` grades each library and ChIP experiment against the ENCODE TF
ChIP-seq audits. The cutoffs, operators, levels, and audit names copy the
portal audit code in [ENCODE-DCC/encoded][encoded] at `9bd8ff1`. Each row
takes a qc.tsv `metric` and `value` and adds a `level` and the `audit` it
raises. Levels rise from `OK` to `WARNING`, `NOT_COMPLIANT`, and `ERROR`.

| Metric | Graded on | Tiers |
|---|---|---|
| `reads_distinct` | GRCh38, hg19, mm10, mm9 | < 5M ERROR, < 10M NOT_COMPLIANT, < 20M WARNING |
| `reads_distinct` | ce10, ce11, dm3, dm6 | < 3M NOT_COMPLIANT, the modERN rule |
| `read_length_mean` | GRCh38, hg19, mm10, mm9 | < 26 ERROR, < 36 NOT_COMPLIANT, < 50 WARNING |
| `nrf` | Every library | < 0.5 NOT_COMPLIANT, < 0.8 WARNING |
| `pbc1` | Every library | < 0.5 NOT_COMPLIANT, < 0.9 WARNING |
| `pbc2` | Every library | < 1 NOT_COMPLIANT, < 10 WARNING. An empty PBC2 with M1 > 0 counts as infinite. |
| `replicates` | ChIP experiments | One replicate is NOT_COMPLIANT. |
| `control` | ChIP experiments | No control is ERROR. |
| `control_reads_distinct` | ChIP experiments with a control, on GRCh38, hg19, mm10, mm9 | The control's summed `reads_distinct`, on the `reads_distinct` tiers |
| `reproducibility` | ChIP experiments with two or more replicates | `borderline` is WARNING. `fail` is NOT_COMPLIANT. |
| `status` | Every experiment | The worst level of the experiment and its libraries, `replicates` and `control` left out |

The first tier whose bound exceeds the value applies.

- Genome keys MUST name their assembly for the count tiers to apply. The
  portal grades counts on no other assembly, so the yeasts get none.
- Usable fragments are `reads_distinct`, the distinct R1 reads on allow-list
  chromosomes. A library without a barcode counts reads as distinct by
  position, as ENCODE does.
- `read_length_mean` is the mean raw read 1 length. The portal grades each
  FASTQ's declared length instead.
- Every target gets the TF tiers. Histone mark classes, NSC, RSC, FRiP,
  antibodies, and metadata go ungraded. ENCODE4 sets no FRiP threshold.
- The sheet has no control type, so any control counts as an input.
- `status` leaves out `replicates` and `control`. A run without controls
  still shows which data pass.

[encoded]: https://github.com/ENCODE-DCC/encoded/tree/9bd8ff18bc5053421d1105a729af4bdb68dab181/src/encoded/audit

## Running

```bash
cp config/config.yaml.template config/config.yaml  # then edit
cp .tests/integration/sra/samples.tsv config/samples.tsv  # then edit
pixi shell                                                 # launcher env
snakemake -n -p --profile profiles/local                   # dry run
snakemake --profile profiles/local                         # local
snakemake --profile profiles/slurm                         # SLURM
```

The first run builds a bowtie2 index per genome under `results/refs/`. GRCh38
took 35 minutes and 6 GB with 16 threads. Later runs reuse it.

Tune `set-resources` in the profiles from the `benchmarks/<run_id>/` TSVs after
a first run.

## Testing

Integration tests follow the Snakemake layout, one case per directory under
`.tests/integration/`. Each case writes its data to a gitignored `data/`.
`.tests/integration/check_run.py` checks any finished run for internal
consistency.

`fixture/make_fixture.py` writes a synthetic run: a random genome, reads around
200 planted sites, and the counts the workflow should report. `check_run.py`
checks a fixture run against those counts too. A second run under another
`run_id` MUST reproduce every tagAlign and narrowPeak, and `--rerun` checks that.

```bash
pixi run --manifest-path workflow/envs/pixi.toml python .tests/integration/fixture/make_fixture.py
snakemake --profile profiles/local --configfile .tests/integration/fixture/data/config.yaml
snakemake --profile profiles/local --configfile .tests/integration/fixture/data/config.yaml \
    --config run_id=fixture_b
pixi run --manifest-path workflow/envs/pixi.toml python .tests/integration/check_run.py \
    results/fixture --expected .tests/integration/fixture/data/expected.json \
    --rerun results/fixture_b
```

`sra/` runs the first 2M reads of 14 SRA runs on dm6, GRCh38, and sacCer3.
They include both layouts, a barcode set, unbarcoded libraries, and a control.
`sra/subsample.py` copies them from the sra-chip store.
`sra/config.yaml` layers over the template, which holds the genomes. The run
needs `/data`. It took 15 minutes on SLURM with the indexes already built.

```bash
pixi run --manifest-path workflow/envs/pixi.toml python .tests/integration/sra/subsample.py
snakemake --profile profiles/slurm \
    --configfile config/config.yaml.template .tests/integration/sra/config.yaml
pixi run --manifest-path workflow/envs/pixi.toml python .tests/integration/check_run.py results/sra
```

`.tests/unit/encode_audit/` holds a QC table with a value on each side of
every audit tier bound, and the audit table it MUST produce.

```bash
pixi run --manifest-path workflow/envs/pixi.toml python workflow/scripts/encode_audit.py \
    --qc .tests/unit/encode_audit/qc.tsv \
    --experiments hs hs_input dm mm sc \
    --genomes GRCh38 GRCh38 dm6 mm10 sacCer3 \
    --controls hs hs_input \
    --output /dev/stdout \
| diff - .tests/unit/encode_audit/encode.tsv
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

- `default`: cutadapt, bowtie2, samtools, bedtools, MACS3, figwig, gawk, and
  polars.
- `idr`: IDR 2.0.4.2, which needs python 3.9 and numpy < 1.20.
- `qc`: matplotlib and poikilos, for the figures.

pixi must be on `PATH`. An `onstart` hook installs the envs before the first
job. Parallel SLURM jobs would otherwise race to install them. After a
manifest edit, run `pixi lock` in `workflow/envs/`. Any change to a manifest
or lock reruns every rule.

## References

Danecek, Petr, James K. Bonfield, Jennifer Liddle, et al. 2021. “Twelve Years of SAMtools and BCFtools.” GigaScience 10 (2): giab008. https://doi.org/10.1093/gigascience/giab008.

He, Qiye, Jeff Johnston, and Julia Zeitlinger. 2015. “ChIP-Nexus Enables Improved Detection of in Vivo Transcription Factor Binding Footprints.” Nature Biotechnology 33 (4): 395–401. https://doi.org/10.1038/nbt.3121.

Landt, Stephen G., Georgi K. Marinov, Anshul Kundaje, et al. 2012. “ChIP-Seq Guidelines and Practices of the ENCODE and modENCODE Consortia.” Genome Research 22 (9): 1813–31. https://doi.org/10.1101/gr.136184.111.

Langmead, Ben, and Steven L. Salzberg. 2012. “Fast Gapped-Read Alignment with Bowtie 2.” Nature Methods 9 (4): 357–59. https://doi.org/10.1038/nmeth.1923.

Li, Qunhua, James B. Brown, Haiyan Huang, and Peter J. Bickel. 2011. “Measuring Reproducibility of High-Throughput Experiments.” The Annals of Applied Statistics 5 (3): 1752–79. https://doi.org/10.1214/11-AOAS466.

Martin, Marcel. 2011. “Cutadapt Removes Adapter Sequences from High-Throughput Sequencing Reads.” EMBnet.Journal 17 (1): 10–12. https://doi.org/10.14806/ej.17.1.200.

Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, et al. 2021. “Sustainable Data Analysis with Snakemake.” F1000Research 10 (January): 33. https://doi.org/10.12688/f1000research.29032.3.

Quinlan, Aaron R., and Ira M. Hall. 2010. “BEDTools: A Flexible Suite of Utilities for Comparing Genomic Features.” Bioinformatics 26 (6): 841–42. https://doi.org/10.1093/bioinformatics/btq033.

Schreiber, Jacob. n.d. Figwig. Github. Accessed October 12, 2026. https://github.com/jmschrei/figwig.

Zhang, Yong, Tao Liu, Clifford A. Meyer, et al. 2008. “Model-Based Analysis of ChIP-Seq (MACS).” Genome Biology 9 (9): R137. https://doi.org/10.1186/gb-2008-9-9-r137.
