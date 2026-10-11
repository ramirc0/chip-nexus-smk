# CLAUDE.md

ENCODE-style ChIP-nexus workflow. README.md holds the steps, sample sheet,
outputs, and style conventions. The behavior choices and their reasons live
in `../COMPARISON.md`, local and untracked. Follow both.

## Commands

```bash
pixi run lint                                # ruff check, ruff format --check, numpydoc lint
pixi run snakemake -n -p --profile profiles/local
pixi run snakemake --profile profiles/slurm
```

## Two pixi workspaces

- Root `pixi.toml` is the launcher: Snakemake and its plugins, plus the `lint` env.
- `workflow/envs/pixi.toml` holds the rule envs: `default` for the tools, `idr` (python 3.9, numpy < 1.20), and `qc` for figures.
- `pyproject.toml` is ruff and numpydoc config only. pixi ignores it.
- Rules use `locked=True`. After a manifest edit, run `pixi lock` in `workflow/envs/`.
- Any manifest or lock change reruns every rule.
- `SOFTWARE_ENV`, `IDR_ENV`, and `QC_ENV` resolve `../envs` relative to the rule file. New rule files MUST sit in `workflow/rules/`.

## Gotchas

- `.smk` files are not linted. Keep helper docstrings numpydoc by hand.
- Every tool a rule calls MUST be in that rule's env. Jobs keep the user's PATH, so a missing tool resolves to a user binary that may need CPU features some nodes lack.
- Cutadapt pass 1 cuts the UMI with `-u`, then matches the fixed barcode anchored at the new start. Pass 2 has no `--rename`, so the pass-1 read ID survives.
- UMI contract: every barcoded read ID ends in `_<KEY>`, the canonical fixed barcode plus the UMI as read. Dedup tools read the UMI from there.
- `dedup.smk` is swappable. A replacement MUST keep the dedup contract below.
- The deliverable BAM is R1 only. bam2bw counts every mapped read with no dup, MAPQ, or mate filter, so it MUST only see deliverable BAMs.
- `chrom_sizes` is the allow-list for tagAligns, tracks, peaks, FRiP, and the dedup histogram. BAMs keep all chromosomes.
- The `rep` wildcard excludes `pooled` with a lookahead. `$` would anchor the whole path, not the wildcard.
- IDR below 20 merged peaks writes empty outputs and status `too_few_peaks`. Other IDR errors fail the job.
- MACS3 caps peaks inside awk. `head` exits early and fails the pipe under pipefail.
- Pseudoreps call `gawk`. mawk gives a different `srand(0)` stream.
- SLURM nodes z1056, z1058, z1062, and z1016 lack `/data`, where the sra-chip FASTQs live. The slurm profile excludes them.
- The launcher's `[activation.env]` sets `OPENBLAS_CORETYPE=Haswell`. That is a setting for this cluster, not workflow logic. Its Broadwell and Cascade Lake nodes otherwise get different OpenBLAS kernels and different MACS3 summits. SLURM jobs inherit it through `--export=ALL`.
- Figures use poikilos `plain` and `pk.save_figure()`, which writes SVG and PNG. Plot rules MUST declare both. Plot scripts read TSVs only. The numbers come from data rules in the `default` env.
- Rules compress and decompress with `python -m zlib_ng.gzip_ng`, never `gzip` or `zcat`. Scripts use `zlib_ng.gzip_ng`. Tools that read `.gz` themselves (bowtie2, MACS3, bedtools, idr, polars) are left alone. It has no `zcat -f`, so FASTQs MUST be gzipped.
- Scripts that import numpy MUST `export OMP_NUM_THREADS={threads}`. Otherwise numpy starts a thread per host core.
- Resources MUST be flat values. Never derive `mem_mb` or `runtime` from input size.
- Auto mode may block `pixi lock`, `pixi install`, or workflow runs, especially for the poikilos git dependency. Ask the user to run them.

## Dedup contract

Inputs:
- `<lib>.filt.bam`, coordinate-sorted and filtered. PE mates carry `fixmate -m` tags.
- The genome's allow-list chrom sizes.
- `dedup_mode()` (`umi`, `position`, or `none`) and the library layout (`se` or `pe`).

Outputs:
- `<lib>.bam` and `.bai`: duplicates removed (kept for `none`), R1 only, all chromosomes.
- `<lib>.dedup.tsv`: `group_size`, `n_groups` over kept R1 reads on allow-list chromosomes. `qc_table.py` derives NRF, PBC1, and PBC2 from it alone.
- Raw tool stats for the user. No rule or script parses them.

## Verifying a change

1. `pixi run lint`.
2. Synthetic fixture, exact expected counts, run twice (commands in README.md "Testing"). `.tests/integration/check_run.py` with `--rerun` MUST report 0 failed.
3. Real subsamples: `.tests/integration/sra/`, the first 2M reads of each run. Commands in README.md "Testing". `check_run.py` MUST report 0 failed.

`resources/`, `results/`, `logs/`, `benchmarks/`, and `.tests/integration/*/data/` are gitignored.
