""" """

from.coordinate import Coordinate
from .sequence import Sequence
from .region import Region
from .variant import VariantRecord, VTYPES
from .utils import adjust_position, compute_indel_length, flatten_list

from typing import Set, List, Dict
from collections import defaultdict

import pandas as pd
import numpy as np

CNAMES = ["haplotype", "variants", "samples"]


class Haplotype(Region):
    def __init__(self, sequence: Sequence, coord: Coordinate, phased: bool, chromcopy: int) -> None:
        super().__init__(sequence, coord)
        self._size = len(self.sequence)
        self._posmap = {i: [i] for i in range(len(self.sequence))}  # positions map
        self._maxpos = self._size - 1
        self._reverse_posmap()  # initialize reverse position map
        self._variants = "NA"
        self._samples = "REF"
        self._samples_suffix = "0/1"  # store the suffix to append to samples 
        if phased:  # if phased, record on which copy the haplotype occurs
            self._samples_suffix = "1|0" if chromcopy == 0 else "0|1"

    def substring(self, start: int, stop: int) -> str:
        start_i, stop_i = self._posmap[start][-1], self._posmap[stop][-1]
        return "".join(self._sequence._sequence_raw[start_i:stop_i])
    
    def _update_map(self, position: int, offset: int, insertion: str) -> None:
        # update position map to handle indel presence in the sequence
        if insertion:  # adjust for insertion
            assert len(self._posmap[position]) == 1
            self._posmap[position] = [self._posmap[position][-1] + i for i in range(offset + 1)]
            for pos_seq, pos_rel in self._posmap.items():
                if pos_seq > position:
                    self._posmap[pos_seq] = [p + offset for p in pos_rel]
        else:  # adjust for deletion
            offset_tracker = 0  # track offset for deletion
            posnew = [self._posmap[position][-1]]  # indel start position
            for pos_seq, pos_rel in self._posmap.items():
                if pos_seq > position:
                    if (offset_tracker := offset_tracker + 1) < offset:
                        self._posmap[pos_seq] = posnew
                    else:
                        self._posmap[pos_seq] = [p - offset for p in pos_rel]

    def _update_sequence(self, start: int, stop: int, alt: str) -> None:
        self._sequence._sequence_raw = self._sequence._sequence_raw[:start] + list(alt.lower()) + self._sequence._sequence_raw[stop:]
    
    def _modify_position(self, position: int, ref: str, alt: str, vtype: str) -> None:
        # compute the variant position within the region sequence
        posrel = adjust_position(position - 1, self._coordinates.start)
        offset = compute_indel_length(ref, alt) if vtype == VTYPES[1] and len(ref) > len(alt) else 0  # compute position offset, if variant is deletion
        # check consistency between reference alleles in sequence and vcf
        if (posrel_stop := posrel + offset + 1) > self._maxpos:
            # avoid out of bounds if indel occurs at boundaries
            posrel_stop = self._maxpos
        idx = posrel_stop - posrel
        sequpdate_offset = len(ref) - len(ref[:idx])
        refnt = self.substring(posrel, posrel_stop)  # retrieve ref nt
        if refnt[:idx] != ref[:idx]:
            raise ValueError(f"Mismatching reference alleles in VCF and reference sequence at position {position} ({refnt} - {ref})")
        # retrieve start and stop positions in relative string
        if any(len(p) > 1 for p in [self._posmap[posrel], self._posmap[posrel_stop]]):
            raise ValueError(f"Forbidden haplotype (position: {position}). Are you modifying a position where an indels occur?")
        start, stop = self._posmap[posrel][-1], self._posmap[posrel_stop][-1] + sequpdate_offset
        self._update_sequence(start, stop, alt)  # update sequence 
        if vtype == VTYPES[1]:  # if indel, update positions map
            offset = compute_indel_length(ref, alt)
            self._update_map(posrel, offset, len(ref) < len(alt))

    def _reverse_posmap(self) -> None:
        # reverse positions map to retrieve reference positions from haplotype 
        # positions
        posmap_rev = defaultdict(list)
        for pos_seq, pos_rel in self._posmap.items():  # populate reversed posmap
            for p in pos_rel:
                posmap_rev[p].append(pos_seq)
        self._posmap_rev = dict(posmap_rev)  # cast to regular dict 


    def add_variants(self, variants: Set[VariantRecord], samples: Set[str]) -> None:
        # sort variants to first add snps, then indels to haplotype
        variants = _sort_variants(variants)
        # modify region sequence to reflect input haplotype
        for variant in variants:  # insert variant in reference sequence
            self._modify_position(variant.position, variant.ref, variant.alt[0], variant.vtype[0])
        # update with current haplotype   
        self._sequence = Sequence("".join(self._sequence._sequence_raw), lower=True)  
        # add variants and samples data associated to current haplotype
        self._variants = ",".join(sorted(flatten_list([v.id for v in variants])))
        self._samples = ",".join(sorted([f"{s}:{self._samples_suffix}" for s in samples]))
        self._reverse_posmap()  # reverse position map

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
    
    @property
    def position_map(self) -> Dict[int, List[int]]:
        return self._posmap
    
    @property
    def position_map_rev(self) -> Dict[int, List[int]]:
        return self._posmap_rev
    
def _sort_variants(variants: Set[VariantRecord]) -> List[VariantRecord]:
    # sort variants set to have snps before indels
    snps, indels = [], []
    for variant in variants:
        if variant.vtype[0] == VTYPES[0]:  # snp
            snps.append(variant)
        else:  # indel
            indels.append(variant)
    return sorted(snps) + sorted(indels)
    

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
        



