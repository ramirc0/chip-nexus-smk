# Run-independent references, built once under results/refs/.


# Every supported genome is under 4 Gbp, so the index is .bt2, never .bt2l.
rule bowtie2_index:
    input:
        lambda w: GENOMES[w.genome]["fasta"],
    output:
        multiext(
            f"{REFDIR}/{{genome}}/{{genome}}",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    params:
        prefix=lambda w: f"{REFDIR}/{w.genome}/{w.genome}",
    log:
        f"{LOGDIR}/bowtie2_index/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/bowtie2_index/{{genome}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bowtie2-build \
            --threads {threads} \
            {input:q} \
            {params.prefix:q}

        sync {output:q}
        """
