""" """

from.coordinate import Coordinate
from .sequence import Sequence
from .region import Region
from .variant import VariantRecord, VTYPES
from .utils import adjust_position, compute_indel_length, flatten_list

from typing import Set, List

import pandas as pd

CNAMES = ["haplotype", "variants", "samples"]


class Haplotype(Region):
    def __init__(self, sequence: Sequence, coord: Coordinate, phased: bool, chromcopy: int) -> None:
        super().__init__(sequence, coord)
        self._size = len(self.sequence)
        self._posmap = {i:i for i in range(len(self.sequence))}  # positions map
        self._samples_suffix = "0/1"  # store the suffix to append to samples 
        if phased:  # if phased, record on which copy the haplotype occurs
            self._samples_suffix = "1|0" if chromcopy == 0 else "0|1"

    def substring(self, start: int, stop: int) -> str:
        start_i, stop_i = self._posmap[start], self._posmap[stop]
        return "".join(self._sequence._sequence_raw[start_i:stop_i])
    
    def _update_map(self, position: int, offset: int) -> None:
        # update position map to handle indel presence in the sequence
        for posorig, posrel in self._posmap.items():
            if posorig > position:
                if (posnew:= posrel + offset) >= self._size - 1:
                    posnew = self._size - 1
                self._posmap[posorig] = posnew


    def _update_sequence(self, start: int, stop: int, alt: str) -> None:
        self._sequence._sequence_raw = self._sequence._sequence_raw[:start] + list(alt) + self._sequence._sequence_raw[stop:]
    
    def _modify_position(self, position: int, ref: str, alt: str, vtype: str) -> None:
        # compute the variant position within the region sequence
        posrel = adjust_position(position - 1, self._coordinates._start)
        print(f"posrel: {posrel} {len(self._sequence._sequence_raw)}")
        print(self._sequence._sequence_raw[posrel:])
        offset = compute_indel_length(ref, alt) if vtype == VTYPES[1] and len(ref) > len(alt) else 0  # compute position offset, if variant is deletion
        # check consistency between reference alleles in sequence and vcf
        refnt = self.substring(posrel, posrel + offset + 1)  # retrieve ref nt
        if refnt != ref:
            raise ValueError(f"Mismatching reference alleles in VCF and reference sequence ({refnt} - {ref})")
        # retrieve start and stop positions in relative string
        start, stop = self._posmap[posrel], self._posmap[posrel + offset + 1]
        if any(isinstance(p, list) for p in [start, stop]):
            raise ValueError(f"Forbidden haplotype (position: {position}). Are you modifying a position where an indels occur?")
        self._update_sequence(start, stop, alt)  # update sequence 
        if vtype == VTYPES[1]:  # if indel, update positions map
            offset = compute_indel_length(ref, alt)
            self._update_map(posrel, offset)


    def add_variants(self, variants: Set[VariantRecord], samples: Set[str]) -> None:
        # sort variants to first add snps, then indels to haplotype
        variants = _sort_variants(variants)
        # modify region sequence to reflect input haplotype
        for variant in variants:  # insert variant in reference sequence
            self._modify_position(variant.position, variant.ref, variant.alt[0], variant.vtype[0])
        # update with current haplotype   
        self._sequence = Sequence("".join(self._sequence._sequence_raw))  
        # add variants and samples data associated to current haplotype
        self._variants = ",".join(sorted(flatten_list([v.id for v in variants])))
        self._samples = ",".join(sorted([f"{s}:{self._samples_suffix}" for s in samples]))

    @property
    def variants(self) -> str:
        if not hasattr(self, "_variants"):
            raise AttributeError(f"The current {self.__class__.__name__} object does not have _variants attribute")
        return self._variants
    
    @property
    def samples(self) -> str:
        if not hasattr(self, "_samples"):
            raise AttributeError(f"The current {self.__class__.__name__} object does not have _samples attribute")
        return self._samples
    
    @property
    def sequence_raw(self) -> List[str]:
        return self._sequence._sequence_raw
    
    @property
    def sequence(self) -> str:
        return self._sequence.sequence
    
def _sort_variants(variants: Set[VariantRecord]) -> List[VariantRecord]:
    # sort variants set to have snps before indels
    snps, indels = [], []
    for variant in variants:
        if variant.vtype[0] == VTYPES[0]:  # snp
            snps.append(variant)
        else:  # indel
            indels.append(variant)
    return snps + indels
    

def haplotypes_table(haplotypes: List[Haplotype]) -> pd.DataFrame:
    if not haplotypes:
        raise ValueError("No input haplotype for table construction")
    # initialize dictionary to use to create the haplotypes table
    df = {cname: [] for cname in CNAMES}
    for haplotype in haplotypes:
        df[CNAMES[0]].append(haplotype.sequence)  # haplotype sequence
        df[CNAMES[1]].append(haplotype.variants)  # variants occurring in the haplotype
        df[CNAMES[2]].append(haplotype.samples)  # samples carrying the haplotype
    return pd.DataFrame(df) 
        



