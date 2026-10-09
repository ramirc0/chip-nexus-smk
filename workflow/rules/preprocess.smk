# Barcode check and removal (cutadapt pass 1), then adapter trimming (pass 2).


# Pass 2 R1 adapters: the paper's adapter cut to the 22 nt prefix that reads
# carry (its BamHI tail is replaced by the TruSeq index primer), plus the
# paper's 1-4 nt 5' truncations. -m 22 -O 4 -e 0.2 are the paper's settings.
rule trim:
    input:
        r1=lambda w: [r1 for r1, _ in LIBRARIES[w.exp, w.rep]["runs"]],
        r2=lambda w: [r2 for _, r2 in LIBRARIES[w.exp, w.rep]["runs"] if r2],
    output:
        fastq=temp(f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.trim.fastq.gz"),
        barcode_json=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.cutadapt.barcode.json",
        adapter_json=f"{OUTDIR}/{{exp}}/{{rep}}/{{exp}}.{{rep}}.cutadapt.adapter.json",
    params:
        layout=lambda w: LIBRARIES[w.exp, w.rep]["layout"],
        barcode=lambda w: barcode_flags(w.exp, w.rep),
        r2_adapters=lambda w: r2_adapter_flags(w.exp, w.rep),
    log:
        f"{LOGDIR}/trim/{{exp}}.{{rep}}.txt",
    benchmark:
        f"{BENCHDIR}/trim/{{exp}}.{{rep}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        if [ {params.layout:q} = pe ]; then
            # The rename puts _<KEY> after the ID, so bowtie2 can no longer
            # strip a /1 /2 suffix and mates would get different names.
            for r in {input.r1:q} {input.r2:q}; do
                id=$(zcat -f "$r" | head -n 1 | cut -d ' ' -f 1 || true)
                case "$id" in
                    */1 | */2)
                        echo "$r: read ID $id ends in /1 or /2" >&2
                        exit 1
                        ;;
                esac
            done
            cutadapt \
                --interleaved \
                -j {threads} \
                {params.barcode:q} \
                --json {output.barcode_json:q} \
                -o - \
                <(zcat -f {input.r1:q}) \
                <(zcat -f {input.r2:q})
        else
            cutadapt \
                -j {threads} \
                {params.barcode:q} \
                --json {output.barcode_json:q} \
                -o - \
                <(zcat -f {input.r1:q})
        fi \
        | cutadapt \
            -j {threads} \
            -a AGATCGGAAGAGCACACGTCTG \
            -a GATCGGAAGAGCACACGTCTG \
            -a ATCGGAAGAGCACACGTCTG \
            -a TCGGAAGAGCACACGTCTG \
            -a CGGAAGAGCACACGTCTG \
            {params.r2_adapters:q} \
            -m 22 \
            -O 4 \
            -e 0.2 \
            -Z \
            --json {output.adapter_json:q} \
            -o {output.fastq:q} \
            -
        """
