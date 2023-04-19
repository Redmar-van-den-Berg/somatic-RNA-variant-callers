include: "common.smk"
include: "vep.smk"


rule all:
    input:
        vardict_tsv=expand(
            "vardict/{sample}.vep.target.tsv", sample=pep.sample_table["sample_name"]
        ),
        varscan_tsv=expand(
            "varscan/{sample}.raw.vcf.gz", sample=pep.sample_table["sample_name"]
        ),


module vep:
    snakefile:
        "vep.smk"
    config:
        config


### VARDICT ###
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


use rule annot from vep as vardict_annotate with:
    input:
        vcf="vardict/{sample}.raw.vcf.gz",
        genome_fasta=config["genome_fasta"],
    output:
        vep="vardict/{sample}.vep.txt",
        stats="vardict/{sample}.vep_stats.txt",


use rule filter_vep from vep as vardict_filter_vep with:
    input:
        vep=rules.vardict_annotate.output.vep,
        ref_id_mapping=config["ref_id_mapping"],
        scr=srcdir("scripts/filter_vep.py"),


use rule vep_table from vep as vardict_vep_table with:
    input:
        vep=rules.vardict_filter_vep.output.filtered,
        scr=srcdir("scripts/vep-table.py"),
    output:
        table="vardict/{sample}.vep.target.tsv",


### varscan ###
rule varscan:
    input:
        bam=get_bam,
        ref=config["genome_fasta"],
        bed=config["bedfile"],
    output:
        vcf="varscan/{sample}.raw.vcf.gz",
    params:
        af=0.05,
    threads: 3
    singularity:
        containers["varscan"]
    shell:
        """
        samtools mpileup \
            --fasta-ref {input.ref} \
            --max-depth 1000000 \
            -s \
            --positions {input.bed} \
            --no-BAQ \
            {input.bam} \
         | grep -vP '\\t\\t' \
         | varscan mpileup2cns \
            --strand-filter 0 \
            --output-vcf 1 \
            --min-var-freq {params.af} \
            --p-value 0.05 \
         | grep -vP '\\t\./\.|\\t0/0' \
         | bgzip -c > {output.vcf}
        """
