# FRiP, the run-level QC table, and its ENCODE audits.


# Reads are their 5' ends, the λ-exonuclease stop.
rule frip:
    input:
        tagalign=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.tagAlign.gz",
        peaks=f"{OUTDIR}/{{exp}}/peaks/{{exp}}.{{rep}}.bfilt.narrowPeak.gz",
        optimal=f"{OUTDIR}/{{exp}}/{{exp}}.optimal.narrowPeak.gz",
    output:
        f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.frip.tsv",
    log:
        f"{LOGDIR}/frip/{{exp}}.{{rep}}.txt",
    benchmark:
        f"{BENCHDIR}/frip/{{exp}}.{{rep}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        ends() {{
            python -m zlib_ng.gzip_ng -cd {input.tagalign:q} \
            | awk 'BEGIN {{ OFS = "\t" }}
                {{ e = ($6 == "+") ? $2 : $3 - 1; print $1, e, e + 1 }}'
        }}

        in_peaks() {{
            ends | bedtools intersect -u -a stdin -b "$1" | wc -l
        }}

        total=$(ends | wc -l)
        {{
            printf 'peak_set\tn_reads\tn_in_peaks\n'
            printf 'replicate\t%s\t%s\n' "$total" "$(in_peaks {input.peaks:q})"
            printf 'optimal\t%s\t%s\n' "$total" "$(in_peaks {input.optimal:q})"
        }} > {output:q}

        sync {output:q}
        """


# Reads per 500 bp bin, the deeptools default. Tracks hold allow-list
# chromosomes only, so bins do too.
rule fingerprint:
    input:
        plus=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.+.bw",
        minus=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.-.bw",
    output:
        f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.fingerprint.tsv",
    log:
        f"{LOGDIR}/fingerprint/{{exp}}.{{rep}}.txt",
    benchmark:
        f"{BENCHDIR}/fingerprint/{{exp}}.{{rep}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # Otherwise numpy starts one thread per core of the host.
        export OMP_NUM_THREADS={threads}

        python workflow/scripts/fingerprint.py \
            --plus {input.plus:q} \
            --minus {input.minus:q} \
            --bin-size 500 \
            --output {output:q}

        sync {output:q}
        """


# Strand profile of each replicate and the pool around the top optimal summits.
rule footprint:
    input:
        peaks=f"{OUTDIR}/{{exp}}/{{exp}}.optimal.narrowPeak.gz",
        plus=lambda w: [f"{lib_prefix(w.exp, s)}.+.bw" for s in samples_of(w.exp)],
        minus=lambda w: [f"{lib_prefix(w.exp, s)}.-.bw" for s in samples_of(w.exp)],
    output:
        f"{OUTDIR}/{{exp}}/{{exp}}.footprint.tsv",
    params:
        samples=lambda w: samples_of(w.exp),
    log:
        f"{LOGDIR}/footprint/{{exp}}.txt",
    benchmark:
        f"{BENCHDIR}/footprint/{{exp}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # Otherwise numpy starts one thread per core of the host.
        export OMP_NUM_THREADS={threads}

        python workflow/scripts/footprint.py \
            --peaks {input.peaks:q} \
            --samples {params.samples:q} \
            --plus {input.plus:q} \
            --minus {input.minus:q} \
            --width 100 \
            --top 1000 \
            --output {output:q}

        sync {output:q}
        """


rule qc_table:
    input:
        barcode=[f"{lib_prefix(e, r)}.cutadapt.barcode.json" for e, r in LIBRARIES],
        adapter=[f"{lib_prefix(e, r)}.cutadapt.adapter.json" for e, r in LIBRARIES],
        bowtie2=[f"{lib_prefix(e, r)}.bowtie2.txt" for e, r in LIBRARIES],
        idxstats=[f"{lib_prefix(e, r)}.filt.idxstats" for e, r in LIBRARIES],
        dedup=[f"{lib_prefix(e, r)}.dedup.tsv" for e, r in LIBRARIES],
        insert_size=[
            f"{lib_prefix(e, r)}.samtools_stats.txt"
            for (e, r), lib in LIBRARIES.items()
            if lib["layout"] == "pe"
        ],
        frip=[f"{lib_prefix(e, r)}.frip.tsv" for e in CHIPS for r in EXPERIMENTS[e]["reps"]],
        peaks=[
            f"{OUTDIR}/{e}/peaks/{e}.{s}.bfilt.narrowPeak.gz" for e in CHIPS for s in samples_of(e)
        ],
        idr_status=[
            f"{OUTDIR}/{e}/idr/{e}.{p}.idr.status" for e in CHIPS for p in idr_pairs(e)
        ],
        idr_peaks=[idr_peaks(e, p) for e in CHIPS for p in idr_pairs(e)],
        reproducibility=[f"{OUTDIR}/{e}/{e}.reproducibility.tsv" for e in CHIPS],
        footprint=[f"{OUTDIR}/{e}/{e}.footprint.tsv" for e in CHIPS],
    output:
        f"{OUTDIR}/qc.tsv",
    log:
        f"{LOGDIR}/qc_table.txt",
    benchmark:
        f"{BENCHDIR}/qc_table.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/qc_table.py \
            --barcode {input.barcode:q} \
            --adapter {input.adapter:q} \
            --bowtie2 {input.bowtie2:q} \
            --idxstats {input.idxstats:q} \
            --dedup {input.dedup:q} \
            --insert {input.insert_size:q} \
            --frip {input.frip:q} \
            --peaks {input.peaks:q} \
            --idr-status {input.idr_status:q} \
            --idr-peaks {input.idr_peaks:q} \
            --reproducibility {input.reproducibility:q} \
            --footprint {input.footprint:q} \
            --output {output:q}

        sync {output:q}
        """


rule encode_audit:
    input:
        f"{OUTDIR}/qc.tsv",
    output:
        f"{OUTDIR}/encode.tsv",
    params:
        experiments=list(EXPERIMENTS),
        genomes=[e["genome"] for e in EXPERIMENTS.values()],
        controls=[x for exp, e in EXPERIMENTS.items() if e["control"] for x in (exp, e["control"])],
    log:
        f"{LOGDIR}/encode_audit.txt",
    benchmark:
        f"{BENCHDIR}/encode_audit.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/encode_audit.py \
            --qc {input:q} \
            --experiments {params.experiments:q} \
            --genomes {params.genomes:q} \
            --controls {params.controls:q} \
            --output {output:q}

        sync {output:q}
        """
