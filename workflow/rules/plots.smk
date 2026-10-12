# QC figures. Plot scripts read only TSVs and QC files. The data rules live in
# qc.smk and peaks.smk. Every figure is an SVG with a PNG beside it.


rule plot_fingerprint:
    input:
        chip=lambda w: [f"{lib_prefix(w.exp, r)}.fingerprint.tsv" for r in EXPERIMENTS[w.exp]["reps"]],
        control=lambda w: control_files(w.exp, "fingerprint.tsv"),
    output:
        svg=f"{OUTDIR}/{{exp}}/{{exp}}.fingerprint.svg",
        png=f"{OUTDIR}/{{exp}}/{{exp}}.fingerprint.png",
    log:
        f"{LOGDIR}/plot_fingerprint/{{exp}}.txt",
    benchmark:
        f"{BENCHDIR}/plot_fingerprint/{{exp}}.tsv"
    software:
        QC_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_fingerprint.py \
            --chip {input.chip:q} \
            --control {input.control:q} \
            --output {output.svg:q}
        """


rule plot_footprint:
    input:
        f"{OUTDIR}/{{exp}}/{{exp}}.footprint.tsv",
    output:
        svg=f"{OUTDIR}/{{exp}}/{{exp}}.footprint.svg",
        png=f"{OUTDIR}/{{exp}}/{{exp}}.footprint.png",
    log:
        f"{LOGDIR}/plot_footprint/{{exp}}.txt",
    benchmark:
        f"{BENCHDIR}/plot_footprint/{{exp}}.tsv"
    software:
        QC_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_footprint.py \
            --input {input:q} \
            --output {output.svg:q}
        """


rule plot_duplicates:
    input:
        lambda w: [f"{lib_prefix(w.exp, r)}.dedup.tsv" for r in EXPERIMENTS[w.exp]["reps"]],
    output:
        svg=f"{OUTDIR}/{{exp}}/{{exp}}.duplicates.svg",
        png=f"{OUTDIR}/{{exp}}/{{exp}}.duplicates.png",
    log:
        f"{LOGDIR}/plot_duplicates/{{exp}}.txt",
    benchmark:
        f"{BENCHDIR}/plot_duplicates/{{exp}}.tsv"
    software:
        QC_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_duplicates.py \
            --input {input:q} \
            --output {output.svg:q}
        """


# Paired-end replicates only.
rule plot_insert_size:
    input:
        lambda w: [f"{lib_prefix(w.exp, r)}.samtools_stats.txt" for r in pe_reps(w.exp)],
    output:
        svg=f"{OUTDIR}/{{exp}}/{{exp}}.insert_size.svg",
        png=f"{OUTDIR}/{{exp}}/{{exp}}.insert_size.png",
    log:
        f"{LOGDIR}/plot_insert_size/{{exp}}.txt",
    benchmark:
        f"{BENCHDIR}/plot_insert_size/{{exp}}.tsv"
    software:
        QC_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_insert_size.py \
            --input {input:q} \
            --output {output.svg:q}
        """


rule plot_idr:
    input:
        f"{OUTDIR}/{{exp}}/idr/{{exp}}.{{pair}}.idr.txt.gz",
    output:
        svg=f"{OUTDIR}/{{exp}}/idr/{{exp}}.{{pair}}.idr.svg",
        png=f"{OUTDIR}/{{exp}}/idr/{{exp}}.{{pair}}.idr.png",
    params:
        threshold=IDR_THRESHOLD,
        title=lambda w: f"{w.exp} {w.pair}",
    log:
        f"{LOGDIR}/plot_idr/{{exp}}.{{pair}}.txt",
    benchmark:
        f"{BENCHDIR}/plot_idr/{{exp}}.{{pair}}.tsv"
    software:
        QC_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_idr.py \
            --input {input:q} \
            --threshold {params.threshold:q} \
            --title {params.title:q} \
            --output {output.svg:q}
        """


rule plot_funnel:
    input:
        f"{OUTDIR}/qc.tsv",
    output:
        svg=f"{OUTDIR}/qc_funnel.svg",
        png=f"{OUTDIR}/qc_funnel.png",
    log:
        f"{LOGDIR}/plot_funnel.txt",
    benchmark:
        f"{BENCHDIR}/plot_funnel.tsv"
    software:
        QC_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_funnel.py \
            --input {input:q} \
            --output {output.svg:q}
        """
