pepfile: config["pepfile"]


containers = {
    "debian": "docker://debian:latest",
    "freebayes": "docker://quay.io/biocontainers/freebayes:1.3.6--h6f59eb7_3",
    "gatk": "docker://quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0",
    "vardict": "docker://quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0",
    "varscan": "docker://quay.io/biocontainers/mulled-v2-58936b48a08c4e06505165c6f560ec9460b431ea:ef260d10ee671f4c7bd8e783939839bb2e0b684e-0",
    "vep": "docker://quay.io/biocontainers/ensembl-vep:108.2--pl5321h4a94de4_0",
    "xopen": "docker://quay.io/biocontainers/xopen:0.7.3--py_0",
}


def get_bam(wildcards):
    """Get the BAM file from the PEP"""
    return pep.sample_table.loc[wildcards.sample, "bamfile"]


def get_vcf(wildcards):
    """Get the VCF file from the PEP"""
    return pep.sample_table.loc[wildcards.sample, "vcffile"]
