# Deduplication. Swap the method by replacing this file. Any replacement MUST
# keep the dedup contract in CLAUDE.md: same inputs, the `_<KEY>` read ID
# suffix as the UMI, and the same deliverable BAM and `.dedup.tsv` outputs.


# `none` keeps duplicates in the BAM, flagged. The histogram counts each
# kept R1 on an allow-list chromosome once, at its `dc` group size.
rule dedup:
    input:
        bam=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.filt.bam",
        sizes=lambda w: genome_file(w.exp, "chrom_sizes"),
    output:
        bam=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.bam",
        bai=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.bam.bai",
        hist=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.dedup.tsv",
        stats=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.markdup.json",
    params:
        barcode=lambda w: (
            ["--barcode-rgx", "_([ACGTN]+)$"] if dedup_mode(w.exp, w.rep) == "umi" else []
        ),
        drop=lambda w: [] if dedup_mode(w.exp, w.rep) == "none" else ["-F", "1024"],
        r1=lambda w: ["-f", "64"] if LIBRARIES[w.exp, w.rep]["layout"] == "pe" else [],
    log:
        f"{LOGDIR}/dedup/{{exp}}.{{rep}}.txt",
    benchmark:
        f"{BENCHDIR}/dedup/{{exp}}.{{rep}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools markdup \
            -@ {threads} \
            --duplicate-count \
            --json \
            -f {output.stats:q} \
            {params.barcode:q} \
            {input.bam:q} \
            - \
        | samtools view \
            -b \
            {params.drop:q} \
            {params.r1:q} \
            --write-index \
            -o {output.bam:q}##idx##{output.bai:q} \
            -

        {{
            printf 'group_size\tn_groups\n'
            samtools view -F 1024 {output.bam:q} \
            | awk 'BEGIN {{ OFS = "\t" }}
                NR == FNR {{ keep[$1]; next }}
                $3 in keep {{
                    size = 1
                    for (i = 12; i <= NF; i++) if ($i ~ /^dc:i:/) size = substr($i, 6)
                    n[size]++
                }}
                END {{ for (s in n) print s, n[s] }}' \
                {input.sizes:q} - \
            | sort -k1,1n
        }} > {output.hist:q}
        """
