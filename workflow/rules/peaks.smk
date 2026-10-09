# Pseudoreplicates, MACS3 peaks, IDR, and the ENCODE conservative and optimal sets.

import math

# ENCODE TF ChIP-seq settings.
PEAK_CAP = 500000
IDR_THRESHOLD = 0.05


# gawk, not awk: mawk gives a different srand(0) stream, so other splits.
rule pseudoreps:
    input:
        f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.tagAlign.gz",
    output:
        pr1=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.pr1.tagAlign.gz",
        pr2=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.pr2.tagAlign.gz",
    log:
        f"{LOGDIR}/pseudoreps/{{exp}}.{{rep}}.txt",
    benchmark:
        f"{BENCHDIR}/pseudoreps/{{exp}}.{{rep}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        workdir=$(mktemp -d)
        trap 'rm -rf "$workdir"' EXIT

        python -m zlib_ng.gzip_ng -cd {input:q} \
        | gawk -v pr1="$workdir/pr1" -v pr2="$workdir/pr2" '
            BEGIN {{ srand(0) }}
            {{ print > (rand() < 0.5 ? pr1 : pr2) }}'
        touch "$workdir/pr1" "$workdir/pr2"
        python -m zlib_ng.gzip_ng -cn "$workdir/pr1" > {output.pr1:q}
        python -m zlib_ng.gzip_ng -cn "$workdir/pr2" > {output.pr2:q}
        """


rule pooled_pseudoreps:
    input:
        pr1=lambda w: [f"{lib_prefix(w.exp, r)}.pr1.tagAlign.gz" for r in EXPERIMENTS[w.exp]["reps"]],
        pr2=lambda w: [f"{lib_prefix(w.exp, r)}.pr2.tagAlign.gz" for r in EXPERIMENTS[w.exp]["reps"]],
    output:
        pr1=f"{OUTDIR}/{{exp}}/pooled/{{exp}}.pooled.pr1.tagAlign.gz",
        pr2=f"{OUTDIR}/{{exp}}/pooled/{{exp}}.pooled.pr2.tagAlign.gz",
    log:
        f"{LOGDIR}/pooled_pseudoreps/{{exp}}.txt",
    benchmark:
        f"{BENCHDIR}/pooled_pseudoreps/{{exp}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        cat {input.pr1:q} > {output.pr1:q}
        cat {input.pr2:q} > {output.pr2:q}
        """


# Sorted by p-value (col 8) and capped inside awk: `head` would exit early
# and fail the pipe. The unfiltered set feeds IDR; `.bfilt` is the deliverable.
rule macs3:
    input:
        treatment=lambda w: (
            f"{OUTDIR}/{w.exp}/{w.peakset.split('.')[0]}/{w.exp}.{w.peakset}.tagAlign.gz"
        ),
        controls=lambda w: control_files(w.exp, "tagAlign.gz"),
        blacklist=lambda w: genome_file(w.exp, "blacklist"),
    output:
        peaks=f"{OUTDIR}/{{exp}}/peaks/{{exp}}.{{peakset}}.narrowPeak.gz",
        bfilt=f"{OUTDIR}/{{exp}}/peaks/{{exp}}.{{peakset}}.bfilt.narrowPeak.gz",
    params:
        control=lambda w: control_flags(w.exp),
        gsize=lambda w: GENOMES[EXPERIMENTS[w.exp]["genome"]]["gsize"],
        shift=lambda w: -round(config["smooth_window"] / 2),
        extsize=lambda w: config["smooth_window"],
        cap=PEAK_CAP,
    log:
        f"{LOGDIR}/macs3/{{exp}}.{{peakset}}.txt",
    benchmark:
        f"{BENCHDIR}/macs3/{{exp}}.{{peakset}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        workdir=$(mktemp -d)
        trap 'rm -rf "$workdir"' EXIT

        macs3 callpeak \
            -t {input.treatment:q} \
            {params.control:q} \
            -f BED \
            -n peaks \
            --outdir "$workdir" \
            -g {params.gsize:q} \
            -p 0.01 \
            --nomodel \
            --shift {params.shift:q} \
            --extsize {params.extsize:q} \
            --keep-dup all \
            --call-summits

        LC_ALL=C sort -k8,8gr "$workdir/peaks_peaks.narrowPeak" \
        | awk -v cap={params.cap:q} 'BEGIN {{ OFS = "\t" }}
            NR <= cap {{ $4 = "Peak_" NR; print }}' \
        | python -m zlib_ng.gzip_ng -cn \
        > {output.peaks:q}

        bedtools intersect -v -a {output.peaks:q} -b {input.blacklist:q} \
        | python -m zlib_ng.gzip_ng -cn \
        > {output.bfilt:q}
        """


# IDR exits non-zero below 20 merged peaks. That case writes empty outputs
# and status `too_few_peaks`, so one weak experiment cannot block qc.tsv.
rule idr:
    input:
        peaks=lambda w: [
            f"{OUTDIR}/{w.exp}/peaks/{w.exp}.{p}.narrowPeak.gz" for p in idr_pairs(w.exp)[w.pair]
        ],
        blacklist=lambda w: genome_file(w.exp, "blacklist"),
    output:
        peaks=f"{OUTDIR}/{{exp}}/idr/{{exp}}.{{pair}}.idr.narrowPeak.gz",
        raw=f"{OUTDIR}/{{exp}}/idr/{{exp}}.{{pair}}.idr.txt.gz",
        status=f"{OUTDIR}/{{exp}}/idr/{{exp}}.{{pair}}.idr.status",
    params:
        threshold=IDR_THRESHOLD,
        min_log10=-math.log10(IDR_THRESHOLD),
    log:
        f"{LOGDIR}/idr/{{exp}}.{{pair}}.txt",
    benchmark:
        f"{BENCHDIR}/idr/{{exp}}.{{pair}}.tsv"
    software:
        IDR_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        workdir=$(mktemp -d)
        trap 'rm -rf "$workdir"' EXIT

        set -- {input.peaks:q}
        if ! idr \
            --samples "$1" "$2" \
            --peak-list "$3" \
            --input-file-type narrowPeak \
            --output-file "$workdir/idr" \
            --rank p.value \
            --soft-idr-threshold {params.threshold:q} \
            --use-best-multisummit-IDR \
            2> "$workdir/stderr"
        then
            cat "$workdir/stderr" >&2
            grep -q 'at least 20 peaks post-merge' "$workdir/stderr"
            : > "$workdir/idr"
            echo too_few_peaks > {output.status:q}
        else
            cat "$workdir/stderr" >&2
            echo ok > {output.status:q}
        fi

        python -m zlib_ng.gzip_ng -cn "$workdir/idr" > {output.raw:q}

        # Column 12 is the global -log10(IDR).
        awk -v min={params.min_log10:q} 'BEGIN {{ OFS = "\t" }}
            $12 >= min {{
                if ($2 < 0) $2 = 0
                print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
            }}' "$workdir/idr" \
        | LC_ALL=C sort -u \
        | LC_ALL=C sort -k8,8gr \
        | bedtools intersect -v -a stdin -b {input.blacklist:q} \
        | python -m zlib_ng.gzip_ng -cn \
        > {output.peaks:q}
        """


rule peak_sets:
    input:
        true_pairs=lambda w: [idr_peaks(w.exp, p) for p in idr_pairs(w.exp) if not p.endswith(".pr")],
        self_pseudo=lambda w: [idr_peaks(w.exp, f"{r}.pr") for r in EXPERIMENTS[w.exp]["reps"]],
        pooled_pseudo=lambda w: (
            [idr_peaks(w.exp, "pooled.pr")] if len(EXPERIMENTS[w.exp]["reps"]) > 1 else []
        ),
    output:
        conservative=f"{OUTDIR}/{{exp}}/{{exp}}.conservative.narrowPeak.gz",
        optimal=f"{OUTDIR}/{{exp}}/{{exp}}.optimal.narrowPeak.gz",
        table=f"{OUTDIR}/{{exp}}/{{exp}}.reproducibility.tsv",
    log:
        f"{LOGDIR}/peak_sets/{{exp}}.txt",
    benchmark:
        f"{BENCHDIR}/peak_sets/{{exp}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/peak_sets.py \
            --true-pairs {input.true_pairs:q} \
            --self-pseudo {input.self_pseudo:q} \
            --pooled-pseudo {input.pooled_pseudo:q} \
            --conservative {output.conservative:q} \
            --optimal {output.optimal:q} \
            --table {output.table:q}
        """
