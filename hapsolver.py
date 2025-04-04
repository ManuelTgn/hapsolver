""" """

from hapsolver import Bed, Fasta, RegionList, VCF, Sequence, HaplotypeGraph, Haplotype, haplotypes_table, flatten_list

from typing import Optional
from time import time

import os


def extract_regions(
    fastafile: str,
    bedfile: str,
    fastaidx: Optional[str] = "",
    padding: Optional[int] = 0,
) -> RegionList:
    """Extract regions from a FASTA file based on a BED file.

    Extracts genomic regions defined in a BED file from a FASTA file, handling
    optional FASTA indexing and padding of regions.

    Args:
        fastafile: Path to the FASTA file.
        bedfile: Path to the BED file.
        fastaidx: Optional path to the FASTA index file.
        padding: Optional padding to add to the start and stop coordinates of each
            region.

    Returns:
        A RegionList object containing the extracted regions.
    """

    fasta = Fasta(fastafile, faidx=fastaidx)  # create fasta object
    bed = Bed(bedfile, padding)  # create bed object to store a set of genomic regions
    # extract regions defined in bed from input fasta
    return bed.extract_regions(fasta)

def read_vcf(regions: RegionList, vcf_file: str):
    vcf = VCF(vcf_file)
    for r in regions:
        variants = vcf.fetch(r.coordinates)
    return variants



if __name__ == "__main__":
    start = time()
    fasta = "/Users/manuel/Desktop/Hapsolver/data/genome/test.fa"
    bed = "/Users/manuel/Desktop/Hapsolver/data/regions/test.bed"
    vcf = "/Users/manuel/Desktop/Hapsolver/data/VCFs/test/chrx.vcf.gz"
    regions = extract_regions(fasta, bed)
    for r in regions:
        print(r)
    variants = read_vcf(regions, vcf)
    vcf = VCF(vcf)
    for v in variants:
        print(v)
    # variants_split = [e for sl in [v.split() for v in variants] for e in sl]
    variants_split = flatten_list([v.split() for v in variants])
    print()
    hg0 = HaplotypeGraph(variants_split)
    hg0.compute_graph()
    # hg0.print()
    print(len(hg0._nodes))
    print()
    haps = hg0.retrieve_haplotypes()
    for h in haps:
        print(h)
    print()
    haps_lst = [Haplotype(Sequence(regions[0].sequence.sequence), regions[0].coordinates, vcf.phased, hg0.chromcopy)]
    for h in haps:
        hap = Haplotype(Sequence(regions[0].sequence.sequence), regions[0].coordinates, vcf.phased, hg0.chromcopy)
        print(h[1])
        hap.add_variants(h[0], h[1])
        haps_lst.append(hap)
    # print(haplotypes_table(haps_lst))
    table = haplotypes_table(haps_lst)
    # for i in table.index:
    #     print(table.iloc[i, 0])
    #     print(table.iloc[i, 1])
    #     print(table.iloc[i, 2])
    #     print()
    print()
    # hg1 = HaplotypeGraph(variants_split, chromcopy=1)
    # hg1.compute_graph()
    # # hg1.print()
    # print(len(hg1._nodes))
    # print()
    # haps = hg1.retrieve_haplotypes()
    # for h in haps:
    #     print(h)
    # print()
    # print()
    print(f"Elapsed time {(time() - start):.6f}s")
