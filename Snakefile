include: "common.smk"
include: "vep.smk"


rule all:
    input:
        vardict_vcf=expand(
            "vardict/{sample}.raw.vcf.gz", sample=pep.sample_table["sample_name"]
        ),
        tsv=expand(
            "vardict/{sample}.vep.target.tsv", sample=pep.sample_table["sample_name"]
        ),


module vep:
    snakefile:
        "vep.smk"
    config:
        config


use rule annot from vep as annotate with:
    input:
        vcf="vardict/{sample}.raw.vcf.gz",
        genome_fasta=config["genome_fasta"],
    output:
        vep="vardict/{sample}.vep.txt",
        stats="vardict/{sample}.vep_stats.txt",


use rule filter_vep from vep as filter_vep with:
    input:
        vep=rules.annotate.output.vep,
        ref_id_mapping=config["ref_id_mapping"],
        scr=srcdir("scripts/filter_vep.py"),


use rule vep_table from vep as vep_table with:
    input:
        vep=rules.filter_vep.output.filtered,
        scr=srcdir("scripts/vep-table.py"),
    output:
        table="vardict/{sample}.vep.target.tsv",


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
