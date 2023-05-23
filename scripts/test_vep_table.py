#!/usr/bin/env python3

import pytest
from vep_table import hgvs_like
from vep_table import info_to_dict, format_to_dict


VARIANTS = [
        # Small insertion
        ("chr1", 114708595, 114708594, "insertion", "-/T", "NC_000001.11:g.114708594_114708595insT"),
        # Big insertion
        ("chr5", 171410540, 171410539, "insertion", "-/TCTG", "NC_000005.10:g.171410539_171410540insTCTG"),
        # SNP
        ("chr1", 114716123, 114716123, "SNV", "C/T", "NC_000001.11:g.114716123C>T"),
        # Small deletion
        ("chr11", 32435164, 32435164, "deletion", "G/-", "NC_000011.10:g.32435164_32435164del"),
        # Big deletion
        ("chr19", 33301849, 33301857, "deletion", "GGCGGCGGC/-", "NC_000019.10:g.33301849_33301857del"),
        # Really big deletion (3748bp), in those cases the allele string is set to "deletion"
        ("chr1", 114705669, 114709417, "deletion", "deletion", "NC_000001.11:g.114705669_114709417del"),
        # Susbstitution
        ("chr4", 54733129, 54733130, "substitution", "GA/TC", "NC_000004.12:g.54733129_54733130delinsTC"),
        # Indel, delete two, insert one base
        ("chr11", 32434944, 32434945, "indel", "GC/T", "NC_000011.10:g.32434944_32434945delinsT"),
        # Indel, delete three, insert one base
        ("chr4", 105236254, 105236256, "indel", "GAG/C", "NC_000004.12:g.105236254_105236256delinsC"),
        # Indel, delete two, insert three bases
        ("chr19", 33302112, 33302113, "indel", "GC/CGT", "NC_000019.10:g.33302112_33302113delinsCGT"),
        # Really big duplication (43504bp), in those cases the allele string is set to "duplicatoin"
        ("chr4", 105190454, 105233958, "duplication", "duplication", "NC_000004.12:g.105190454_105233958dup"),
]

NOT_IMPLEMENTED_VARIANTS = [
    # An insertion where the reference allele is not '-'
    ("chr1", 114708595, 114708594, "insertion", "T/-", ""),
    # A SNV where start is not equal to end
    ("chr1", 114716123, 1, "SNV", "C/T", "NC_000001.11:g.114716123C>T"),
    # A deletion where the alt allele is not '-'
    ("chr11", 32435164, 32435164, "deletion", "-/G", ""),
    # A substitution where ref and alt are of different size
    ("chr1", 1, 2, "substitution", "G/TC", ""),
    # Unsupported variant classes
    ("chr1", 1, 2, "sequence_alteration", "G/TC", ""),
    ("chr1", 1, 2, "sequence_variation", "G/TC", ""),
]


INFO = [
    ("SAMPLE=sample1;DP=5", {"SAMPLE": "sample1", "DP": "5"}),
]

FORMAT = [
    ("GT:DP", "0/1:5", {"GT": "0/1", "DP": "5"}),
]


@pytest.mark.parametrize(["chrom", "start", "end", "variant_class", "allele_string", "hgvs"], VARIANTS)
def test_hgvs_like(chrom, start, end, variant_class, allele_string, hgvs):
    assert hgvs_like(chrom, start, end, variant_class, allele_string) == hgvs


@pytest.mark.parametrize(["chrom", "start", "end", "variant_class", "allele_string", "hgvs"], NOT_IMPLEMENTED_VARIANTS)
def test_no_implementedhgvs_like(chrom, start, end, variant_class, allele_string, hgvs):
    with pytest.raises(NotImplementedError):
        hgvs_like(chrom, start, end, variant_class, allele_string)


@pytest.mark.parametrize(["info", "result"], INFO)
def test_info_field(info, result):
    assert info_to_dict(info) == result


@pytest.mark.parametrize(["header", "fields", "result"], FORMAT)
def test_info_field(header, fields, result):
    assert format_to_dict(header, fields) == result
