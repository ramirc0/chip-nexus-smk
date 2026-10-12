# Raw 5' count tracks per strand, both positive.


# figwig bam2bw counts every mapped read: no duplicate, MAPQ, or mate filter.
# The BAMs are already filtered and hold R1 only. It drops chromosomes absent
# from -s.
rule tracks:
    input:
        bams=lambda w: [f"{lib_prefix(w.exp, r)}.bam" for r in reps_of(w.exp, w.sample)],
        sizes=lambda w: genome_file(w.exp, "chrom_sizes"),
    output:
        plus=f"{OUTDIR}/{{exp}}/{{sample}}/{{exp}}.{{sample}}.+.bw",
        minus=f"{OUTDIR}/{{exp}}/{{sample}}/{{exp}}.{{sample}}.-.bw",
    params:
        name=lambda w: lib_prefix(w.exp, w.sample),
    log:
        f"{LOGDIR}/tracks/{{exp}}.{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/tracks/{{exp}}.{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        figwig bam2bw \
            -s {input.sizes:q} \
            -n {params.name:q} \
            -p {threads} \
            {input.bams:q}
        """
