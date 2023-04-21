include: "common.smk"
include: "vep.smk"


rule all:
    input:
        vardict_tsv=expand(
            "{caller}/{sample}.vep.target.tsv", sample=pep.sample_table["sample_name"], caller=config["callers"]
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
    output:
        filtered="vardict/{sample}.vep.target.txt.gz",


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
    container:
        containers["varscan"]
    log:
        "log/varscan.{sample}.txt",
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
            --p-value 0.05 2> {log} \
         | grep -vP '\\t\./\.|\\t0/0' \
         | bgzip -c > {output.vcf}
        """


use rule annot from vep as varscan_annotate with:
    input:
        vcf="varscan/{sample}.raw.vcf.gz",
        genome_fasta=config["genome_fasta"],
    output:
        vep="varscan/{sample}.vep.txt",
        stats="varscan/{sample}.vep_stats.txt",


use rule filter_vep from vep as varscan_filter_vep with:
    input:
        vep=rules.varscan_annotate.output.vep,
        ref_id_mapping=config["ref_id_mapping"],
        scr=srcdir("scripts/filter_vep.py"),
    output:
        filtered="varscan/{sample}.vep.target.txt.gz",


use rule vep_table from vep as varscan_vep_table with:
    input:
        vep=rules.varscan_filter_vep.output.filtered,
        scr=srcdir("scripts/vep-table.py"),
    output:
        table="varscan/{sample}.vep.target.tsv",


### Mutect2 ###
rule split_N_CIGAR:
    input:
        bam=get_bam,
        ref=config["genome_fasta"],
        ref_dict=config["genome_dict"],
        bed=config["bedfile"],
    output:
        bam="mutect2/{sample}.split.bam",
        bai="mutect2/{sample}.split.bai",
    log:
        "log/split_CIGAR.{sample}.txt",
    threads: 8
    container:
        containers["gatk"]
    shell:
        """
        gatk SplitNCigarReads \
            --input {input.bam} \
            --reference {input.ref} \
            --output {output.bam} \
            --intervals {input.bed} \
            --create-output-bam-index 2> {log}
        """


rule mutect2:
    input:
        bam=rules.split_N_CIGAR.output.bam,
        bai=rules.split_N_CIGAR.output.bai,
        ref=config["genome_fasta"],
        pon=config["pon"],
    output:
        vcf="mutect2/{sample}.raw.vcf.gz",
    threads: 8
    log:
        "log/mutect2.{sample}.txt",
    container:
        containers["gatk"]
    shell:
        """
        gatk Mutect2 \
            --reference {input.ref} \
            --input {input.bam} \
            --panel-of-normals {input.pon} \
            --output {output.vcf} 2> {log}
        """


use rule annot from vep as mutect2_annotate with:
    input:
        vcf="mutect2/{sample}.raw.vcf.gz",
        genome_fasta=config["genome_fasta"],
    output:
        vep="mutect2/{sample}.vep.txt",
        stats="mutect2/{sample}.vep_stats.txt",


use rule filter_vep from vep as mutect2_filter_vep with:
    input:
        vep=rules.mutect2_annotate.output.vep,
        ref_id_mapping=config["ref_id_mapping"],
        scr=srcdir("scripts/filter_vep.py"),
    output:
        filtered="mutect2/{sample}.vep.target.txt.gz",


use rule vep_table from vep as mutect2_vep_table with:
    input:
        vep=rules.mutect2_filter_vep.output.filtered,
        scr=srcdir("scripts/vep-table.py"),
    output:
        table="mutect2/{sample}.vep.target.tsv",


### Freebayes ###
rule freebayes:
    input:
        bam=get_bam,
        ref=config["genome_fasta"],
        bed=config["bedfile"],
    params:
        af=0.05,
    output:
        vcf="freebayes/{sample}.raw.vcf.gz",
    threads: 8
    log:
        "log/freebayes.{sample}.txt",
    container:
        containers["freebayes"]
    shell:
        """
        freebayes \
            --fasta-reference {input.ref} \
            --bam {input.bam} \
            --min-alternate-fraction {params.af} \
            --targets {input.bed} \
            --pooled-continuous \
            2> {log} | bgzip > {output.vcf}
        """


use rule annot from vep as freebayes_annotate with:
    input:
        vcf="freebayes/{sample}.raw.vcf.gz",
        genome_fasta=config["genome_fasta"],
    output:
        vep="freebayes/{sample}.vep.txt",
        stats="freebayes/{sample}.vep_stats.txt",


use rule filter_vep from vep as freebayes_filter_vep with:
    input:
        vep=rules.freebayes_annotate.output.vep,
        ref_id_mapping=config["ref_id_mapping"],
        scr=srcdir("scripts/filter_vep.py"),
    output:
        filtered="freebayes/{sample}.vep.target.txt.gz",


use rule vep_table from vep as freebayes_vep_table with:
    input:
        vep=rules.freebayes_filter_vep.output.filtered,
        scr=srcdir("scripts/vep-table.py"),
    output:
        table="freebayes/{sample}.vep.target.tsv",
