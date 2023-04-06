include: "common.smk"


rule all:
    input:
        vardict_vcf=expand(
            "vardict/{sample}.raw.vcf.gz", sample=pep.sample_table["sample_name"]
        ),


rule vardict:
    input:
        bam=get_bam,
        ref=config["genome_fasta"],
        bed=config["bedfile"],
    output:
        vcf="vardict/{sample}.raw.vcf.gz",
    params:
        af=0.05,
        bed_chrom=1,
        bed_start=2,
        bed_end=3,
        bed_gene=4,
    log:
        "log/vardict.{sample}.txt",
    threads: 3
    container:
        containers["vardict"]
    shell:
        """
        vardict-java \
            -G {input.ref} \
            -N {wildcards.sample} \
            -b {input.bam} \
            -f {params.af} \
            -c {params.bed_chrom} \
            -S {params.bed_start} \
            -E {params.bed_end} \
            -g {params.bed_gene} \
            --verbose \
            {input.bed} \
            |
        teststrandbias.R \
            |
        var2vcf_valid.pl -A | gzip > {output.vcf} 2> {log}
        """
