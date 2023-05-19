include: "common.smk"


rule all_vep:
    input:
        vep=expand("{sample}/vep.target.tsv", sample=pep.sample_table["sample_name"]),


rule annot:
    """Annotate variants using VEP

    This rule is not used directly, but intended to be modified and reused in other pipelines
    """
    input:
        vcf="placeholder",
        genome_fasta="placeholder",
    params:
        online="" if config.get("cache_vep") else "--database",
        offline=f" --offline --cache_version 108 --everything --merged --dir {config['cache_vep']}"
        if config.get("cache_vep")
        else "",
        freq_filter=" --af_gnomade --check_frequency" if config.get("cache_vep") else "",
        max_pop_af=config.get("max_pop_af", 0.01),
    output:
        vep="{sample}.vep.placeholder",
        stats="{sample}.stats.placeholder",
    log:
        "log/annotate.{sample}.txt",
    threads: 8
    container:
        containers["vep"]
    shell:
        """
        vep \
            -i {input.vcf} \
            --fasta {input.genome_fasta} \
            {params.online} \
            {params.offline} \
            {params.freq_filter} \
            --fork {threads} \
            --allele_number --stats_text --json --force_overwrite --assembly GRCh38 \
            --format vcf \
            --polyphen b \
            --sift b \
            --hgvs \
            --af \
            --freq_pop gnomADe \
            --freq_freq {params.max_pop_af} \
            --freq_gt_lt gt \
            --freq_filter exclude \
            --stats_file {output.stats} \
            --output_file STDOUT | gzip > {output.vep}\
            2> {log}
        """


use rule annot as annotate with:
    input:
        vcf=get_vcf,
        genome_fasta=config["genome_fasta"],
    output:
        vep="{sample}/vep.txt.gz",
        stats="{sample}/vep_stats.txt",


rule filter_vep:
    input:
        vep=rules.annotate.output.vep,
        ref_id_mapping=config["ref_id_mapping"],
        scr=srcdir("scripts/filter_vep.py"),
    params:
        vep_consequences=config["vep_include_consequences"],
    output:
        filtered="{sample}/vep.target.txt.gz",
    log:
        "log/filter_vep.{sample}.txt",
    threads: 1
    container:
        containers["vep"]
    shell:
        """
        python {input.scr} \
            {input.vep} \
            {input.ref_id_mapping} \
            --consequences {params.vep_consequences} \
            | gzip > {output.filtered} 2> {log}
        """


rule vep_table:
    input:
        vep=rules.filter_vep.output.filtered,
        scr=srcdir("scripts/vep_table.py"),
        cache=config.get("hgvs_cache", list()),
    params:
        cache=f"--hgvs-cache-in {config['hgvs_cache']}"
        if config.get("hgvs_cache")
        else "",
    output:
        table="{sample}/vep.target.tsv",
        cache="{sample}/hgvs.cache",
    log:
        "log/vep_table.{sample}.txt",
    container:
        containers["xopen"]
    shell:
        """
        python {input.scr} \
            {params.cache} \
            --hgvs-cache-out {output.cache} \
            {input.vep} > {output.table} 2> {log}
        """
