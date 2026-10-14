# Alignment and read filtering, then per-library and pooled tagAligns.


# MAPQ >= 30 and -F 1804 follow ENCODE. bowtie2 writes mates next to each
# other, so fixmate needs no name sort. It clears the proper-pair flag of
# mates whose partner was filtered, and the second -f 2 drops them.
# --reorder keeps the order of tied reads, and so the pseudoreps, independent
# of thread timing.
rule align:
    input:
        fastq=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.trim.fastq.gz",
        index=lambda w: multiext(
            bowtie2_index(w.exp),
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    output:
        bam=temp(f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.filt.bam"),
        idxstats=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.filt.idxstats",
        bowtie2=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.bowtie2.txt",
    params:
        index=lambda w: bowtie2_index(w.exp),
        layout=lambda w: LIBRARIES[w.exp, w.rep]["layout"],
        reads=lambda w: (
            ["-X", "2000", "--dovetail", "--interleaved"]
            if LIBRARIES[w.exp, w.rep]["layout"] == "pe"
            else ["-U"]
        ),
        proper=lambda w: ["-f", "2"] if LIBRARIES[w.exp, w.rep]["layout"] == "pe" else [],
    log:
        f"{LOGDIR}/align/{{exp}}.{{rep}}.txt",
    benchmark:
        f"{BENCHDIR}/align/{{exp}}.{{rep}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        workdir=$(mktemp -d)
        trap 'rm -rf "$workdir"' EXIT

        if [ {params.layout:q} = pe ]; then
            drop_orphans() {{
                samtools fixmate -m -u - - | samtools view -u -f 2 -
            }}
        else
            drop_orphans() {{
                cat
            }}
        fi

        bowtie2 \
            --end-to-end \
            --reorder \
            -p {threads} \
            -x {params.index:q} \
            {params.reads:q} {input.fastq:q} \
            2> {output.bowtie2:q} \
        | samtools view -u -q 30 -F 1804 {params.proper:q} - \
        | drop_orphans \
        | samtools sort \
            -@ {threads} \
            -m 1G \
            -T "$workdir/sort" \
            -o {output.bam:q} \
            -

        # R1 per chromosome, the reference count for dedup and track checks.
        samtools view -u -F 128 {output.bam:q} \
        | samtools idxstats - \
        > {output.idxstats:q}

        sync {output:q}
        """


# Insert sizes of the filtered pairs, before dedup keeps R1 only.
rule insert_size:
    input:
        f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.filt.bam",
    output:
        f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.samtools_stats.txt",
    log:
        f"{LOGDIR}/insert_size/{{exp}}.{{rep}}.txt",
    benchmark:
        f"{BENCHDIR}/insert_size/{{exp}}.{{rep}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools stats {input:q} > {output:q}

        sync {output:q}
        """


# One rule for replicates and pools: the input is one BAM, or every replicate
# BAM for `pooled`. Keeps allow-list chromosomes only. PE BAMs hold R1 only.
rule tagalign:
    input:
        bams=lambda w: [f"{lib_prefix(w.exp, r)}.bam" for r in reps_of(w.exp, w.sample)],
        sizes=lambda w: genome_file(w.exp, "chrom_sizes"),
    output:
        f"{OUTDIR}/{{exp}}/{{sample}}/{{exp}}.{{sample}}.tagAlign.gz",
    log:
        f"{LOGDIR}/tagalign/{{exp}}.{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/tagalign/{{exp}}.{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        for bam in {input.bams:q}; do
            bedtools bamtobed -i "$bam"
        done \
        | awk 'BEGIN {{ OFS = "\t" }}
            NR == FNR {{ keep[$1]; next }}
            $1 in keep {{ print $1, $2, $3, "N", 1000, $6 }}' \
            {input.sizes:q} - \
        | python -m zlib_ng.gzip_ng -cn \
        > {output:q}

        sync {output:q}
        """
