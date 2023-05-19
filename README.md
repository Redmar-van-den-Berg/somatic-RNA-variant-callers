[![Continuous Integration](https://github.com/Redmar-van-den-Berg/somatic-RNA-variant-callers/actions/workflows/ci.yml/badge.svg)](https://github.com/Redmar-van-den-Berg/somatic-RNA-variant-callers/actions/workflows/ci.yml)
[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pep.databio.org)
![GitHub release](https://img.shields.io/github/v/release/redmar-van-den-berg/somatic-RNA-variant-callers)
![Commits since latest release](https://img.shields.io/github/commits-since/redmar-van-den-berg/somatic-RNA-variant-callers/latest)

# somatic-RNA-variant-callers
Compare different somatic variant caller to determine how they perform on
RNAseq data.

## Background
The goal is to find a variant caller that can:
- Call somatic variants
- Does not rely on a matched normal sample
- Can use on RNAseq data

## Installation
Download the repository from github
```bash
git clone https://github.com/Redmar-van-den-Berg/somatic-RNA-variant-callers.git
```

Install and activate the
[conda](https://docs.conda.io/en/latest/miniconda.html)
environment.
```bash
conda env create --file environment.yml
conda activate somatic-RNA-variant-callers
```

## Settings
You can pass settings to the pipeline using the `--config` flag, or in the
specified `--configfile`.

### Supported settings
The following settings are available for the pipeline.
| Option               | Type              | Explanation                             |
| ---------------------| ----------------- | --------------------------------------- |
| genome_fasta         | Required file     | Reference file to use with variant calling |
| genome_fai           | Required file     | .fai index file for the reference file  |
| genome_dict          | Required file     | .dict index file for the reference      |
| pon                  | Required file     | Panel of normals for mutect2            |
| pon_tbi              | Required file     | Index file for panel of normals         |
| bedfile              | Required file     | Bed file with regions to call variants  |
| ref_id_mapping       | Required file     | File with gene/transcripts of interest  |
| callers              | Required list     | One or more of `vardict`, `varscan`, `mutect2`, `freebayes` |
| vep_include_consequences | Required list | [VEP consequences](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html) to include |
| min_af               | Optional float    | Minimum sample allele frequency to call ( 0.05 ) |
| max_pop_af           | Optional float    | Maximum population allele frequency ( 0.01 ) |
| hgvs_cache           | Optional file     | File with cached HGVS normalisations |

## Tests
You can run the tests that accompany this pipeline with the following commands

```bash
# Check if requirements are installed, and run linting on the Snakefile
pytest --kwd --tag sanity

# Test the pipeline settings in dry-run mode
pytest --kwd --tag dry-run

# Test the performance of the pipeline by running on the test data
pytest --kwd --tag integration
```
