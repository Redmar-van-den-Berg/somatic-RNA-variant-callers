include: "common.smk"
include: "vep.smk"


rule all:
    input:
        final_tsv=expand(
            "{caller}/{sample}.vep.target.tsv",
            sample=pep.sample_table["sample_name"],
            caller=config["callers"],
        ),
        final=expand("{caller}/merged.tsv", caller=config["callers"]),


module vep:
    snakefile:
        "vep.smk"
    config:
        config


rule tmpdir:
    output:
        directory("tmp"),
    container:
        containers["varscan"]
    log:
        "log/tmpdir.txt",
    shell:
        """
        mkdir -p {output} 2> {log}
        """


rule merge_table:
    """Merge final VEP tables

    This rule is not used directly, but intended to be modified for each caller
    """
    input:
        tables=list(),
        src=srcdir("scripts/merge_table.py"),
    params:
        samples=get_samples(),
    output:
        "",
    container:
        containers["cutadapt"]
    log:
        "",
    shell:
        """
        python3 {input.src} \
            --samples {params.samples} \
            --tables {input.tables} > {output} 2> {log}
        """


### VARDICT ###
rule vardict:
    input:
        bam=get_bam,
        ref=config["genome_fasta"],
        bed=config["bedfile"],
        tmp=rules.tmpdir.output,
    output:
        vcf="vardict/{sample}.raw.vcf.gz",
    params:
        min_af=config["min_af"],
        bed_chrom=1,
        bed_start=2,
        bed_end=3,
        bed_gene=4,
    log:
        vardict="log/vardict.{sample}.txt",
        strandbias="log/strandbias.{sample}.txt",
        var2vcf="log/var2vcf.{sample}.txt",
        gzip="log/gzip.{sample}.txt",
    threads: 3
    container:
        containers["vardict"]
    shell:
        """
        export TMPDIR={input.tmp}

        vardict-java \
            -G {input.ref} \
            -N {wildcards.sample} \
            -b {input.bam} \
            -f {params.min_af} \
            -c {params.bed_chrom} \
            -S {params.bed_start} \
            -E {params.bed_end} \
            -g {params.bed_gene} \
            --verbose \
            {input.bed} 2> {log.vardict} \
            |
        teststrandbias.R 2> {log.strandbias} \
            |
        var2vcf_valid.pl -A 2> {log.var2vcf} \
            |
        gzip > {output.vcf} 2> {log.gzip}
        """


use rule annot from vep as vardict_annotate with:
    input:
        vcf="vardict/{sample}.raw.vcf.gz",
        genome_fasta=config["genome_fasta"],
    output:
        vep="vardict/{sample}.vep.txt.gz",
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
        scr=srcdir("scripts/vep_table.py"),
    output:
        table="vardict/{sample}.vep.target.tsv",
        cache="vardict/{sample}.hgvs.cache",


use rule merge_table as vardict_merge_table with:
    input:
        tables=[f"vardict/{sample}.vep.target.tsv" for sample in get_samples()],
        src=srcdir("scripts/merge_table.py"),
    output:
        "vardict/merged.tsv",
    log:
        "log/vardict_merge_table.txt",


### varscan ###
rule varscan:
    input:
        bam=get_bam,
        ref=config["genome_fasta"],
        bed=config["bedfile"],
    output:
        vcf="varscan/{sample}.raw.vcf.gz",
    params:
        min_af=config["min_af"],
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
            --min-var-freq {params.min_af} \
            --p-value 0.05 2> {log} \
         | grep -vP '\\t\./\.|\\t0/0' \
         | bgzip -c > {output.vcf}
        """


use rule annot from vep as varscan_annotate with:
    input:
        vcf="varscan/{sample}.raw.vcf.gz",
        genome_fasta=config["genome_fasta"],
    output:
        vep="varscan/{sample}.vep.txt.gz",
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
        scr=srcdir("scripts/vep_table.py"),
    output:
        table="varscan/{sample}.vep.target.tsv",
        cache="varscan/{sample}.hgvs.cache",


use rule merge_table as varscan_merge_table with:
    input:
        tables=[f"varscan/{sample}.vep.target.tsv" for sample in get_samples()],
        src=srcdir("scripts/merge_table.py"),
    output:
        "varscan/merged.tsv",
    log:
        "log/varscan_merge_table.txt",


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
    params:
        min_af=config["min_af"],
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
            --minimum-allele-fraction {params.min_af} \
            --input {input.bam} \
            --panel-of-normals {input.pon} \
            --output {output.vcf} 2> {log}
        """


use rule annot from vep as mutect2_annotate with:
    input:
        vcf="mutect2/{sample}.raw.vcf.gz",
        genome_fasta=config["genome_fasta"],
    output:
        vep="mutect2/{sample}.vep.txt.gz",
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
        scr=srcdir("scripts/vep_table.py"),
    output:
        table="mutect2/{sample}.vep.target.tsv",
        cache="mutect2/{sample}.hgvs.cache",


use rule merge_table as mutect2_merge_table with:
    input:
        tables=[f"mutect2/{sample}.vep.target.tsv" for sample in get_samples()],
        src=srcdir("scripts/merge_table.py"),
    output:
        "mutect2/merged.tsv",
    log:
        "log/mutect2_merge_table.txt",


### Freebayes ###
rule freebayes:
    input:
        bam=get_bam,
        ref=config["genome_fasta"],
        bed=config["bedfile"],
    params:
        min_af=config["min_af"],
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
            --min-alternate-fraction {params.min_af} \
            --targets {input.bed} \
            --pooled-continuous \
            2> {log} | bgzip > {output.vcf}
        """


use rule annot from vep as freebayes_annotate with:
    input:
        vcf="freebayes/{sample}.raw.vcf.gz",
        genome_fasta=config["genome_fasta"],
    output:
        vep="freebayes/{sample}.vep.txt.gz",
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
        scr=srcdir("scripts/vep_table.py"),
    output:
        table="freebayes/{sample}.vep.target.tsv",
        cache="freebayes/{sample}.hgvs.cache",


use rule merge_table as freebayes_merge_table with:
    input:
        tables=[f"freebayes/{sample}.vep.target.tsv" for sample in get_samples()],
        src=srcdir("scripts/merge_table.py"),
    output:
        "freebayes/merged.tsv",
    log:
        "log/freebayes_merge_table.txt",
