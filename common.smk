pepfile: config["pepfile"]


containers = {
    "debian": "docker://debian:latest",
    "vardict": "docker://quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0",
}


def get_bam(wildcards):
    """Get the BAM file from the PEP"""
    return pep.sample_table.loc[wildcards.sample, "bam"]
