""" """

from .coordinate import Coordinate

from typing import Optional, List, Tuple, Set
from pysam import TabixFile, tabix_index

import warnings
import sys
import os


TBI = "tbi"  # tabix index suffix
VTYPES = ["snp", "indel"]  # variant types


class VariantRecord:
    """Represents a variant record from a VCF file.

    Stores information about a variant, including its chromosome, position,
    reference and alternative alleles, variant type, and sample genotypes.
    """

    def __init__(self) -> None:
        """Initialize a VariantRecord object.

        Initializes an empty VariantRecord object. Attributes for the variant
        record will be populated later using the read_vcf_line() method.
        """
        pass

    def __repr__(self) -> str:
        """Return a string representation of the variant record.

        Returns a string representation of the VariantRecord object, including the
        class name, chromosome, position, reference allele, and alternative alleles.

        Returns:
            The string representation of the variant record.
        """
        altalleles = ",".join(self._alt)
        return f'<{self.__class__.__name__} object; variant="{self._chrom} {self._position} {self._ref} {altalleles}">'

    def __str__(self) -> str:
        """Return a string representation of the variant record.

        Returns a string representation of the VariantRecord object, suitable for
        writing to a VCF file. The string includes the chromosome, position,
        reference allele, and alternative alleles, separated by tabs.

        Returns:
            The string representation of the variant record.
        """
        altalleles = ",".join(self._alt)
        return f"{self._chrom}\t{self._position}\t{self._ref}\t{altalleles}"

    def __eq__(self, vrecord: "VariantRecord") -> bool:
        """Check if two VariantRecord objects are equal.

        Two VariantRecord objects are considered equal if they have the same
        chromosome, position, reference allele, and alternative alleles.

        Args:
            vrecord: The other VariantRecord object to compare to.

        Returns:
            True if the two VariantRecord objects are equal, False otherwise.

        Raises:
            AttributeError: If the comparison fails due to missing attributes.
        """
        if not hasattr(vrecord, "_chrom"):
            raise AttributeError(
                f"Comparison between {self.__class__.__name__} object failed"
            )
        if not hasattr(vrecord, "_position"):
            raise AttributeError(
                f"Comparison between {self.__class__.__name__} object failed"
            )
        if not hasattr(vrecord, "_ref"):
            raise AttributeError(
                f"Comparison between {self.__class__.__name__} object failed"
            )
        if not hasattr(vrecord, "_alt"):
            raise AttributeError(
                f"Comparison between {self.__class__.__name__} object failed"
            )
        return (
            self._chrom == vrecord.contig
            and self._position == vrecord.position
            and self._ref == vrecord.ref
            and self._alt == vrecord.alt
        )
    
    def __lt__(self, vrecord: "VariantRecord") -> bool:
        return self._position < vrecord.position
    
    def __gt__(self, vrecord: "VariantRecord") -> bool:
        return self._position > vrecord.position

    def __hash__(self) -> int:
        return hash((self._chrom, self._position, self._ref, tuple(self._alt)))

    def _retrieve_alt_alleles(self, altalleles: str) -> List[str]:
        """Retrieve and parse alternative alleles from a string.

        Retrieves alternative alleles from a comma-separated string.

        Args:
            altalleles: A string containing comma-separated alternative alleles.

        Returns:
            A list of alternative alleles.

        Raises:
            AttributeError: If the input is not a string.
            Exception: For any other unexpected error during parsing.
        """
        if not altalleles:
            return []
        # alternative alleles in multiallelic sites are separated by a comma
        try:
            return altalleles.split(",")
        except AttributeError as e:
            raise AttributeError(
                "Alternative alleles must be encoded in a string"
            ) from e
        except Exception as e:
            raise Exception("Unexpected error while parsing alternative alleles") from e

    def _assess_vtype(self) -> List[str]:
        """Determine the variant type for all alternative alleles.

        Determines if each alternative allele represents a SNP or an indel based on the
        lengths of the reference and alternative alleles.

        Returns:
            A list of variant types (either "snp" or "indel") corresponding to each
            alternative allele.
        """
        assert hasattr(self, "_ref")
        assert hasattr(self, "_alt")
        return [_assign_vtype(self._ref, altallele) for altallele in self._alt]

    def _assign_id(self) -> List[str]:
        """Assign or compute variant IDs.

        Assigns or computes variant IDs based on whether the variant is monoallelic or
        multiallelic. For monoallelic variants, it uses the provided ID or computes one
        if not available. For multiallelic variants, it computes an ID for each
        alternative allele.

        Args:
            vid: The variant ID (may be empty for monoallelic variants).

        Returns:
            A list of variant IDs.
        """
        if self._allelesnum == 1:
            # variant id not available, construct the id using chrom, position, ref,
            # and alt (e.g. chrx-100-A/G)
            return [_compute_id(self._chrom, self._position, self._ref, self._alt[0])]
        # if multiallelic site compute the id for each alternative allele
        # avoid potential confusion due to alternative alleles at same position
        # labeled with the same id
        return [
            _compute_id(self._chrom, self._position, self._ref, altallele)
            for altallele in self._alt
        ]

    def _copy(self, i: int) -> "VariantRecord":
        """Create a copy of the variant record for a specific allele.

        Creates a copy of the current VariantRecord object, representing a single
        alternative allele specified by the index `i`.  This is useful for
        handling multiallelic sites where each alternative allele needs to be
        treated as a separate variant.

        Args:
            i: The index of the alternative allele to copy.

        Returns:
            A new VariantRecord object representing the specified alternative allele.
        """
        # copy current variant record instance
        vrecord = VariantRecord()  # create new instance
        # adjust ref/alt alleles and positions for multiallelic sites
        ref, alt, position = _adjust_multiallelic(
            self._ref, self._alt[i], self._position
        )
        vrecord._chrom = self._chrom
        vrecord._position = position
        vrecord._ref = ref
        vrecord._alt = [alt]
        vrecord._allelesnum = 1
        vrecord._vtype = [self._vtype[i]]
        vrecord._filter = self._filter
        vrecord._vid = [self._vid[i]]
        vrecord._samples = [self._samples[i]]
        return vrecord

    def read_vcf_line(
        self, variant: List[str], samples: List[str], phased: bool
    ) -> None:
        """Read and parse a VCF line.

        Parses a line from a VCF file and populates the attributes of the
        VariantRecord object with the extracted information.

        Args:
            variant: A list of strings representing the fields of the VCF line.
            samples: A list of sample names.
            phased: True if the genotypes are phased, False otherwise.
        """

        self._chrom = variant[0]  # store chromosome
        self._position = int(variant[1])  # store variant position
        self._ref = variant[3]  # store ref allele
        self._alt = self._retrieve_alt_alleles(variant[4])  # store alt alleles
        self._allelesnum = len(self._alt)  # number of alt alleles
        self._vtype = self._assess_vtype()  # establish whether is a snp or indel
        self._filter = variant[6]  # store filter value
        self._vid = self._assign_id()  # assign variant id
        self._samples = _genotypes_to_samples(
            variant[9:], samples, self._allelesnum, phased
        )  # recover samples with their genotypes

    def split(self, vtype: Optional[str] = None) -> List["VariantRecord"]:
        """Split a multiallelic variant record by variant type.

        Splits a multiallelic VariantRecord object into a list of VariantRecord objects,
        each representing a single alternative allele of the specified variant type.

        Args:
            vtype: The variant type to select ("snp" or "indel"). If None, all
                variant types are included.

        Returns:
            A list of VariantRecord objects, one for each alternative allele matching
            the specified variant type.
        """
        vtypes_filter = VTYPES if vtype is None else [vtype]
        return [
            self._copy(i)
            for i, _ in enumerate(self._vtype)
            if self._vtype[i] in vtypes_filter
        ]

    def get_altalleles(self, vtype: str) -> List[str]:
        """Retrieve alternative alleles of a specific variant type.

        Returns a list of alternative alleles that match the specified variant type
        (either "snp" or "indel").

        Args:
            vtype: The variant type to select ("snp" or "indel").

        Returns:
            A list of alternative alleles matching the specified type.
        """
        assert vtype in VTYPES
        # return the alternative alleles representing snps or indels
        return [
            altallele
            for i, altallele in enumerate(self._alt)
            if self._vtype[i] == vtype
        ]

    @property
    def filter(self) -> str:
        return self._filter

    @property
    def contig(self) -> str:
        return self._chrom

    @property
    def position(self) -> int:
        return self._position

    @property
    def ref(self) -> str:
        return self._ref

    @property
    def alt(self) -> List[str]:
        return self._alt

    @property
    def vtype(self) -> List[str]:
        return self._vtype

    @property
    def samples(self) -> Tuple[Set[str], Set[str]]:
        return self._samples

    @property
    def id(self) -> List[str]:
        return self._vid

    @property
    def allelesnum(self) -> int:
        return self._allelesnum


def _assign_vtype(ref: str, alt: str) -> bool:
    """Determine the variant type.

    Determines if a variant is an indel or a SNP based on the lengths of the
    reference and alternative alleles.

    Args:
        ref: The reference allele.
        alt: The alternative allele.

    Returns:
        "indel" if the lengths of the reference and alternative alleles are not equal,
        "snp" otherwise.
    """
    return VTYPES[1] if len(ref) != len(alt) else VTYPES[0]


def _compute_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Compute a variant ID.

    Computes a variant ID using the chromosome, position, reference allele, and
    alternative allele, following the IGVF consortium notation.

    Args:
        chrom: The chromosome.
        pos: The position.
        ref: The reference allele.
        alt: The alternative allele.

    Returns:
        The computed variant ID.
    """
    # compute variant id for variants without id, or multiallelic sites
    # use IGVF consortium notation
    return f"{chrom}-{pos}-{ref}/{alt}"


def _adjust_multiallelic(ref: str, alt: str, pos: int) -> Tuple[str, str, int]:
    """Adjust reference/alternative alleles and position for multiallelic sites.

    Adjusts the reference and alternative alleles, and the variant position for
    multiallelic sites based on the lengths of the original reference and
    alternative alleles.  This function helps normalize variant representation
    for easier comparison and processing.

    Args:
        ref: The original reference allele.
        alt: The original alternative allele.
        pos: The original variant position.

    Returns:
        A tuple containing the adjusted reference allele, alternative allele, and
        variant position.
    """

    if len(ref) == len(alt):  # likely snp
        ref_new, alt_new = ref[-1], alt[-1]  # adjust ref/alt alleles
        pos_new = pos + len(ref) - 1  # ref/alt have same length
    elif len(ref) > len(alt):  # deletion
        ref_new = ref[len(alt) - 1 :]  # adjust ref allele
        alt_new = alt[-1]  # adjust alt allele
        pos_new = pos + (len(alt)) - 1  # adjust variant position
    else:  # insertion
        ref_new = ref[-1]  # adjust ref allele
        alt_new = alt[len(ref) - 1 :]  # adjust alt allele
        pos_new = pos + len(ref) - 1  # adjust variant position
    return ref_new, alt_new, pos_new


def _genotypes_to_samples(
    genotypes: List[str], samples: List[str], allelesnum: int, phased: bool
) -> Tuple[Set[str], Set[str]]:
    """Extract sample information from genotypes.

    Parses genotype strings to determine which samples carry each alternative allele.
    Handles both phased and unphased genotypes.

    Args:
        genotypes: A list of genotype strings.
        samples: A list of sample names.
        allelesnum: The number of alternative alleles.
        phased: True if the genotypes are phased, False otherwise.

    Returns:
        A tuple containing two lists of sets. The first list contains sets of samples
        with the variant on the left copy (or the only copy if unphased), and the
        second list contains sets of samples with the variant on the right copy (only
        relevant for phased data).

    Raises:
        TypeError: If the genotype string cannot be split.
        Exception: If an unexpected error occurs during genotype parsing.
    """
    # define two sets storing samples with variant occurrence on left and right
    # copy respectively
    # if unphased vcf, is used only the left set
    sampleshap = [(set(), set()) for _ in range(allelesnum)]
    gtsep = "|" if phased else "/"  # define genotype separator char
    for i, gt in enumerate(genotypes):
        try:
            gt1, gt2 = gt.split(":")[0].split(gtsep)
        except TypeError as e:
            raise TypeError(f"Split object is not of {str.__name__} type") from e
        except Exception as e:
            raise Exception(
                f"An unexpected error occurred while parsing genotype {gt}"
            ) from e
        if gt1 not in ["0", "."]:  # left copy
            sampleshap[int(gt1) - 1][0].add(samples[i])
        if phased and gt2 not in ["0", "."]:  # right copy
            sampleshap[int(gt2) - 1][1].add(samples[i])
        elif gt2 not in ["0", "."]:  # unphased, so use left copy
            sampleshap[int(gt2) - 1][0].add(samples[i])
    return sampleshap


class VCF:
    def __init__(self, fname: str, vcfidx: Optional[str] = "") -> None:
        if not os.path.isfile(fname):
            raise FileNotFoundError(f"Connot find input VCF {fname}")
        self._fname = fname  # store input filename
        self._vcfidx = self._search_index(vcfidx)  # initialize vcf index
        # initialize TabixFile object with the previously computed index
        self._vcf = TabixFile(self._fname, index=self._vcfidx)
        if len(set(self._vcf.contigs)) != 1:  # assume vcf data about one contig
            raise ValueError(
                f"Input VCF {fname} store variants belonging to multiple contigs"
            )
        self._contig = self._vcf.contigs[0]
        self._samples = self._vcf.header[-1].strip().split()[9:]  # recover samples
        self._phased = False  # by default treat vcf as unphased
        self._is_phased()  # check if the input vcf is phased

    def _search_index(self, vcfidx: Optional[str] = "") -> str:
        """Search for or validate a Tabix index for the VCF file.

        Searches for a Tabix index (.tbi) for the associated VCF file if one is not
        provided. If a path to an index is provided, it validates that the index
        exists and is not empty.

        Args:
            vcfidx: An optional path to a Tabix index file.

        Returns:
            The path to the Tabix index file, or an empty string if not found.

        Raises:
            FileNotFoundError: If the provided index file does not exist or is empty.
        """

        # look for index for the current vcf, if not found compute it
        if not vcfidx:
            if _find_tbi(self._fname):  # index found, store it
                return f"{self._fname}.{TBI}"
            # index not found -> compute it de novo and store it in the same folder
            # as the input vcf
            sys.stdout.write(f"Tabix index not found for {self._fname}\n")
            return ""
        # precomputed vcf index index must be a non empty file
        if not (os.path.isfile(vcfidx) and os.stat(vcfidx).st_size > 0):
            raise FileNotFoundError(f"Not existing or empty VCF index {vcfidx}")
        return vcfidx

    def index_vcf(self, pytest: Optional[bool] = False) -> None:
        """Create or update the Tabix index for the VCF file.

        Creates or updates the Tabix index (.tbi) for the associated VCF file.
        If an index already exists, it will be overwritten.

        Raises:
            OSError: If an error occurs during indexing.
            RuntimeWarning: If an index already exists.
        """
        if self._vcfidx and not pytest:  # launch warning
            warnings.warn("Tabix index already present, forcing update", RuntimeWarning)
        # compute tabix index if not provided during object initialization
        try:  # create index in the same folder as the input vcf
            tabix_index(self._fname, preset="vcf", force=True)
        except OSError as e:
            raise OSError(f"An error occurred while indexing {self._fname}") from e
        assert _find_tbi(self._fname)
        self._faidx = f"{self._fname}.{TBI}"

    def _is_phased(self) -> None:
        """Check if the VCF file is phased.

        Checks if the VCF file is phased by examining the genotype of the first variant.
        If the genotype contains a pipe character '|', the VCF is considered phased.
        The result is stored in the _phased attribute.
        """
        assert hasattr(self, "_vcf")  # otherwise we couldn't establish phasing
        for variant in self._vcf.fetch():  # fecth only the first variant
            gt = variant.strip().split()[9]
            break  # no further iterations required
        # establish from genotype whther the vcf is phased or not
        if "|" in gt:
            self._phased = True

    def fetch(self, coordinate: Coordinate) -> List[VariantRecord]:
        """Fetch variants within a specified genomic interval.

        Fetches variants from the VCF file that overlap with the specified genomic
        interval. Handles potential mismatches between chromosome prefixes and checks for
        phasing.

        Args:
            coordinate: A Coordinate object specifying the genomic interval.

        Returns:
            A list of VariantRecord objects representing the variants within the
            specified interval.

        Raises:
            ValueError: If the VCF and coordinate contigs mismatch or if an invalid
                reference or position is provided.
            IndexError: If an attempt is made to fetch data outside the bounds of the
                indexed file.
            Exception: For any other unexpected error during fetching.
        """

        if self._contig != coordinate.contig:
            # may be just caused by a prefix in contig
            vcfcontig = self._contig.replace("chr", "")
            coordcontig = coordinate.contig.replace("chr", "")
            if vcfcontig != coordcontig:
                raise ValueError(
                    f"Mismatching VCF and coordinate contigs ({self._contig} - {coordinate.contig})"
                )
        try:  # extract variants in the input range from vcf file
            self._is_phased()  # assess whether the vcf is phased
            return [
                _create_variant_record(v.strip().split(), self._samples, self._phased)
                for v in self._vcf.fetch(
                    self._contig, coordinate.start, coordinate.stop
                )
            ]
        except ValueError as e:
            raise ValueError(
                f"Invalid reference or position provided ({self._contig}\t{coordinate.start}\t{coordinate.stop})"
            ) from e
        except IndexError as e:
            raise IndexError(
                f"Tried to fetch data outside the bounds of the indexed file ({self._contig}\t{coordinate.start}\t{coordinate.stop})"
            ) from e
        except Exception as e:
            raise Exception(
                f"An unexpected error occurred ({self._contig}\t{coordinate.start}\t{coordinate.stop})"
            ) from e

    @property
    def contig(self) -> str:
        return self._contig if self._contig.startswith("chr") else f"chr{self._contig}"

    @property
    def phased(self) -> bool:
        return self._phased


def _find_tbi(vcf: str) -> bool:
    """Check if a Tabix index exists for a VCF file.

    Checks if a Tabix index (.tbi) exists for the given VCF file and is a non-empty
    file.

    Args:
        vcf: The path to the VCF file.

    Returns:
        True if the index exists and is a non-empty file, False otherwise.
    """
    # avoid unexpected crashes due to file location
    vcfindex = f"{os.path.abspath(vcf)}.{TBI}"
    if os.path.exists(vcfindex):  # index must be a non empty file
        return os.path.isfile(vcfindex) and os.stat(vcfindex).st_size > 0
    return False


def _create_variant_record(
    variant: List[str], samples: List[str], phased: bool
) -> VariantRecord:
    vrecord = VariantRecord()  # create variant record instance
    vrecord.read_vcf_line(variant, samples, phased)  # read vcf line
    return vrecord
