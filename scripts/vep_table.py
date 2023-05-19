#!/usr/bin/env python3

import argparse
import xopen
import json
import urllib.request as request


FIELDS = [
    "seq_region_name",
    "start",
    "end",
    "variant_class",
    "allele_string",
    "hgvsg"
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
    'chrM': 'NC_012920.1',
    # Stupid online version of VEP renames chrM to MT
    'MT': 'NC_012920.1'
}


def hgvs_like(chrom, start, end, variant_class, allele_string):
    """Create a non-normalised HGVS description for a given variant"""
    chrom = MAPPING[chrom]

    # Variant representation in case of error
    var = f"{variant_class} {chrom}:{start} {allele_string}"

    # Pull out the ref and alt from the allele_string
    try:
        ref, alt = allele_string.split("/")
    except ValueError as e:
        # Big deletions have the word "deletion" as allele string
        if not allele_string == "deletion":
            raise e
        alt = "-"
        ref = None


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


def mutalyzer_normalise(hgvs_like):
    """Use Mutalyzer to normalise the HGVS description"""
    url=f"https://mutalyzer.nl/api/normalize/{hgvs_like}"

    with request.urlopen(url) as response:
        js = json.loads(response.read())
    return js["normalized_description"]
    try:
        with urllib.request.urlopen(url) as response:
            js = json.loads(response.read())
    except:
        print(f"Error for {url}")
        return None
    return js["normalized_description"]


def make_hgvs(variant, cache):
    chrom = variant["seq_region_name"]
    start = variant["start"]
    end = variant["end"]
    variant_class = variant["variant_class"]
    allele_string = variant["allele_string"]

    hgvsg_like = hgvs_like(chrom, start, end, variant_class, allele_string)

    # Add to the cache if it is missing
    if hgvsg_like not in cache:
        hgvsg = mutalyzer_normalise(hgvsg_like)
        cache[hgvsg_like] = hgvsg

    # Always return from the cache
    return cache[hgvsg_like]


def read_vep(fname):
    with xopen.xopen(fname) as fin:
        for line in fin:
            yield json.loads(line)


def read_cache(cache_in):
    """Read the HGVS normalisation cache file"""
    cache = dict()
    if not cache_in:
        return cache

    with open(cache_in) as fin:
        for line in fin:
            key, value = line.strip().split()
            cache[key] = value

    return cache


def write_cache(cache, fname):
    """Write the HGVS normalisation cache file"""
    if not fname:
        return

    with open(fname, "wt") as fout:
        for key, value in cache.items():
            print(key, value, file=fout)


def main(vepfile, sep, cache_in, cache_out):
    # Read in the HGVS cache
    cache = read_cache(cache_in)

    # Print the header
    print(*FIELDS, *TRANSCRIPT_FIELDS, sep=sep)

    # For every variant, print every transcript consequence
    for variant in read_vep(vepfile):
        # Determine HGVS if variant_class is defined (only for
        # offline mode)
        if "variant_class" in variant:
            variant["hgvsg"] = make_hgvs(variant, cache)
        for transcript in variant["transcript_consequences"]:
            # Get all the data we want to print in a list
            variant_data = [variant.get(f) for f in FIELDS]
            transcript_data = [transcript.get(f) for f in TRANSCRIPT_FIELDS]
            print(*variant_data, *transcript_data, sep=sep)

    # Write out the HGVS cache
    write_cache(cache, cache_out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('vep', help='VEP output file')
    parser.add_argument('--sep', default='\t')
    parser.add_argument('--hgvs-cache-in', default=None)
    parser.add_argument('--hgvs-cache-out', default=None)

    args = parser.parse_args()

    if args.hgvs_cache_in and args.hgvs_cache_in == args.hgvs_cache_out:
        msg = f"The same cache ({args.hgvs_cache_in}) cannot be used as input and output"
        raise RuntimeError(msg)

    main(args.vep, args.sep, args.hgvs_cache_in, args.hgvs_cache_out)
