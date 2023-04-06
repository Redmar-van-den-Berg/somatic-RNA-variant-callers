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
| bedfile              | Required file     | Bed file with regions to call variants  |

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
