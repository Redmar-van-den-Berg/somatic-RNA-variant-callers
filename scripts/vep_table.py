#!/usr/bin/env python

import argparse
import xopen
import json


FIELDS = [
    "seq_region_name",
    "start",
    "end",
    "variant_class",
    "allele_string"
]

TRANSCRIPT_FIELDS = [
    "variant_allele",
    "gene_id",
    "gene_symbol",
    "hgvsc",
    "hgvsp",
    "impact",
    "consequence_terms"
]

MAPPING = {
    'chr1': 'NC_000001.11',
    'chr2': 'NC_000002.12',
    'chr3': 'NC_000003.12',
    'chr4': 'NC_000004.12',
    'chr5': 'NC_000005.10',
    'chr6': 'NC_000006.12',
    'chr7': 'NC_000007.14',
    'chr8': 'NC_000008.11',
    'chr9': 'NC_000009.12',
    'chr10': 'NC_000010.11',
    'chr11': 'NC_000011.10',
    'chr12': 'NC_000012.12',
    'chr13': 'NC_000013.11',
    'chr14': 'NC_000014.9',
    'chr15': 'NC_000015.10',
    'chr16': 'NC_000016.10',
    'chr17': 'NC_000017.11',
    'chr18': 'NC_000018.10',
    'chr19': 'NC_000019.10',
    'chr20': 'NC_000020.11',
    'chr21': 'NC_000021.9',
    'chr22': 'NC_000022.11',
    'chrX': 'NC_000023.11',
    'chrY': 'NC_000024.10',
    'chrM': 'NC_012920.1'
}


def hgvs_like(chrom, start, end, variant_class, allele_string):
    """Create a non-normalised HGVS description for a given variant"""
    chrom = MAPPING[chrom]

    # Variant representation in case of error
    var = f"{variant_class} {chrom}:{start} {allele_string}"

    # Pull out the ref and alt from the allele_string
    ref, alt = allele_string.split("/")

    if variant_class == "insertion":
        if not ref == "-":
            raise NotImplementedError(var)
        return f"{chrom}:g.{end}_{start}ins{alt}"

    if variant_class == "SNV":
        if start != end:
            raise NotImplementedError(var)
        return f"{chrom}:g.{start}{ref}>{alt}"

    if variant_class == "deletion":
        if not alt == "-":
            raise NotImplementedError(var)
        return f"{chrom}:g.{start}_{end}del"

    if variant_class == "substitution":
        if len(ref) != len(alt):
            raise NotImplementedError(var)
        return f"{chrom}:g.{start}_{end}delins{alt}"

    if variant_class == "indel":
        return f"{chrom}:g.{start}_{end}delins{alt}"

    # Unhandle variant_class
    raise NotImplementedError(var)


def read_vep(fname):
    with xopen.xopen(fname) as fin:
        for line in fin:
            yield json.loads(line)


def main(fnames, sep):
    print(*FIELDS, *TRANSCRIPT_FIELDS, sep=sep)
    for vep in fnames:
        for variant in read_vep(vep):
            for transcript in variant["transcript_consequences"]:
                # Get all the data we want to print in a list
                variant_data = [variant.get(f) for f in FIELDS]
                transcript_data = [transcript.get(f) for f in TRANSCRIPT_FIELDS]
                print(*variant_data, *transcript_data, sep=sep)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('vep', nargs='+')
    parser.add_argument('--sep', default='\t')

    args = parser.parse_args()

    main(args.vep, args.sep)
