pepfile: config["pepfile"]


containers = {
    "debian": "docker://debian:latest",
    "vardict": "docker://quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0",
    "vep": "docker://quay.io/biocontainers/ensembl-vep:108.2--pl5321h4a94de4_0",
    "xopen": "docker://quay.io/biocontainers/xopen:0.7.3--py_0",
}


def get_bam(wildcards):
    """Get the BAM file from the PEP"""
    return pep.sample_table.loc[wildcards.sample, "bamfile"]


def get_vcf(wildcards):
    """Get the VCF file from the PEP"""
    return pep.sample_table.loc[wildcards.sample, "vcffile"]
