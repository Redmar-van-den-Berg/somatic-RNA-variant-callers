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
    "hgvsc",
    "hgvsp",
    "impact",
    "consequence_terms"
]


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
