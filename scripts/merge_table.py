#!/usr/bin/env python3

import argparse

HEADER = ['seq_region_name', 'start', 'end', 'variant_class', 'allele_string',
          'hgvsg', 'variant_allele', 'gene_id', 'gene_symbol', 'hgvsc',
          'hgvsp', 'impact', 'consequence_terms']

def main(samples, tables):
    # First, we print the header
    print("sample", *HEADER, sep="\t")

    # Next, we prepend the sample name, and print each row
    for s, t in zip(samples, tables):
        with open(t) as fin:
            # Check and skip header
            assert next(fin).strip('\n').split('\t') == HEADER
            for line in fin:
                print(s, line.strip('\n'), sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples', nargs='+')
    parser.add_argument('--tables', nargs='+')

    args = parser.parse_args()

    if len(args.samples) != len(args.tables):
        raise RuntimeError("Unequal number of samples and tables")

    main(args.samples, args.tables)
