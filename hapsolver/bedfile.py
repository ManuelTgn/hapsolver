"""Provides classes and functions for working with BED files.

This module defines the Bed class for representing genomic regions from a BED file,
along with associated iterator and parsing functions. It allows for reading,
manipulating, and extracting regions from BED files, and integrating with FASTA
files for sequence retrieval.
"""


from coordinate import Coordinate
from sequence import Fasta
from region import Region, RegionList

from typing import List, Union

import os


class Bed:
    """Represents genomic regions from a BED file.

    Stores genomic coordinates from a BED file and provides methods for accessing,
    manipulating, and extracting regions.

    Attributes:
        _fname (str): The path to the BED file.
        _coordinates (List[Coordinate]): A list of Coordinate objects representing the genomic regions.
    """

    def __init__(self, bedfile: str, padding: int) -> None:
        """Initialize a Bed object.

        Initializes a Bed object from a BED file, storing genomic coordinates.
        The coordinates are padded by the specified amount.

        Args:
            bedfile: Path to the BED file.
            padding: The amount of padding to add to the start and stop coordinates.

        Raises:
            FileNotFoundError: If the BED file is not found.
        """
        if not os.path.isfile(bedfile):
            raise FileNotFoundError(f"Cannot find input BED file {bedfile}")
        self._fname = bedfile  # store input file name
        # read input bed file content and store a list of coordinates
        self._coordinates = self._read(padding)

    def __repr__(self) -> str:
        """Return a string representation of the Bed object.

        The representation includes the class name and the number of stored regions.

        Returns:
            A string representation of the Bed object.
        """
        return f"<{self.__class__.__name__} object; stored regions={len(self)}>"

    def __len__(self) -> int:
        """Return the number of regions in the Bed object.

        Returns the number of regions stored in the _coordinates attribute.

        Returns:
            The number of regions.

        Raises:
            AttributeError: If the _coordinates attribute is missing.
        """
        if not hasattr(self, "_coordinates"):  # always trace this error
            raise AttributeError(
                f"Missing _coordinates attribute on {self.__class__.__name__}"
            )
        return len(self._coordinates)

    def __iter__(self) -> "BedIterator":
        """Return an iterator over the regions in the Bed object.

        Returns an iterator that allows iterating over the regions in the Bed object.

        Returns:
            An iterator over the regions in the Bed object.
        """
        return BedIterator(self)

    def __getitem__(self, idx: Union[int, slice]) -> Coordinate:
        """Return the regions coordinate at the given index or slice.

        Returns the regions coordinate at the given index or slice from the
        _coordinates list.

        Args:
            idx: The index or slice of the desired regions coordinate.

        Returns:
            The regions coordinate at the given index or slice.

        Raises:
            AttributeError: If the _coordinates attribute is missing.
            IndexError: If the index is out of range.
            TypeError: If the index is not an integer or a slice.
        """
        if not hasattr(self, "_coordinates"):  # always trace this error
            raise AttributeError(
                f"Missing _ccordinates attribute on {self.__class__.__name__}"
            )
        try:  # access _coordinates list to return the corresponding coordinate
            return self._coordinates[idx]
        except IndexError as e:
            raise IndexError(f"Index {idx} out of range") from e
        except TypeError as e:
            raise TypeError(
                f"Invalid index type ({type(idx).__name__}), expected {int.__name__} or {slice.__name__}"
            ) from e

    def _read(self, padding: int) -> List[Coordinate]:
        """Read the BED file and return a list of Coordinate objects.

        Reads the BED file specified by the filename attribute, parses each line,
        and returns a list of Coordinate objects representing the genomic intervals.

        Args:
            padding: The padding to add to the start and stop coordinates.

        Returns:
            A list of Coordinate objects.

        Raises:
            FileNotFoundError: If the BED file is not found.
            PermissionError: If there is a permission error reading the BED file.
            IOError: If there is an I/O error reading the BED file.
            Exception: If any other error occurs while reading the BED file.
        """
        coordinates = []  # list of coordinates from the input bed (Coordinate objs)
        try:
            with open(self._fname, mode="r") as infile:  # begin Bedfile parsing
                # valid bed files must have at least three columns: chromosome, start,
                # and end coordinates; separator is not necessarily tab
                coordinates.extend(
                    _parse_bed_line(line, i + 1, padding)
                    for i, line in enumerate(infile)
                    if not line.startswith("#") and line.strip()
                )
        except FileNotFoundError as e:  # bed file not found
            raise FileNotFoundError(f"Unable to find {self._fname}") from e
        except PermissionError as e:  # permission error on reading
            raise PermissionError(
                f"Permission denied when trying reading {self._fname}"
            ) from e
        except IOError as e:  # i/o error on read
            raise IOError(f"I/O error while reading {self._fname}") from e
        except Exception as e:  # generic exception caught
            raise Exception(
                f"An unexpected error occurred while reading {self._fname}"
            ) from e  # sourcery skip: raise-specific-error
        return coordinates

    def extract_regions(self, fasta: Fasta) -> RegionList:
        """Extract and return a list of regions from a FASTA file.

        Extracts regions from the provided FASTA file based on the coordinates
        stored in the Bed object.

        Args:
            fasta: The Fasta object representing the FASTA file.

        Returns:
            A RegionList object containing the extracted regions.

        Raises:
            AttributeError: If the _coordinates attribute is missing or None.
        """
        if not hasattr(self, "_coordinates") or self._coordinates is None:
            raise AttributeError("Missing coordinates list, cannot extract sequences")
        return RegionList([Region(fasta.fetch(c), c) for c in self._coordinates])


class BedIterator:
    """Iterator for Bed objects.

    Iterates over the regions in a Bed object.

    Attributes:
        bed: The Bed object to iterate over.
        _index: The current index of the iterator.
    """

    def __init__(self, bed: Bed) -> None:
        """Initialize the BedIterator.

        Initializes the iterator with the given Bed object.

        Args:
            bed: The bed object to iterate over.

        Raises:
            AttributeError: If the bed object is missing the _coordinates attribute.
        """
        if not hasattr(bed, "_coordinates"):  # always trace this error
            raise AttributeError(
                f"Missing _ccordinates attribute on {bed.__class__.__name__}"
            )
        self._bed = bed  # bed object to be iterated
        self._index = 0  # iterator index used over the coordinates list

    def __next__(self) -> Coordinate:
        """Return the next region coordinate in the Bed object.

        Returns the next region coordinate in the Bed object, and advances the
        iterator.

        Returns:
            The next region coordinate in the Bed object.

        Raises:
            StopIteration: If there are no more regions to iterate over.
        """
        if self._index < len(self._bed):  # not reached end of bed coordinates
            result = self._bed[self._index]
            self._index += 1  # go to next position in the list
            return result
        raise StopIteration  # stop iteration over bed object


def _parse_bed_line(bedline: str, linenum: int, padding: int) -> Coordinate:
    """Parse a line from a BED file and return a Coordinate object.

    Parses a line from a BED file, extracts the chromosome, start, and stop
    coordinates, and returns a Coordinate object. The start and stop coordinates
    are padded by the given padding value.

    Args:
        bedline: The line from the BED file.
        linenum: The line number in the BED file.
        padding: The padding to add to the start and stop coordinates.

    Returns:
        A Coordinate object representing the genomic interval.

    Raises:
        ValueError: If the line has less than three columns or if the stop coordinate
            is less than the start coordinate.
        TypeError: If the start or stop values are not integers.
    """
    columns = bedline.strip().split()  # recover bed fields for current line
    # minimum fields required: chrom, start, stop
    # (see https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
    if len(columns) < 3:
        raise ValueError(f"Less than three columns at line {linenum}")
    try:  # initialize chrom, start and stop fields
        chrom, start, stop = columns[0], int(columns[1]), int(columns[2])
    except ValueError as e:  # raise if start or stop are not valid int
        raise TypeError(
            f"Start/stop values at line {linenum} are not {int.__name__}"
        ) from e
    if stop < start:  # ensure valid genomic range
        raise ValueError(
            f"Stop < start coordinate ({stop} < {start}) at line {linenum}"
        )
    # if required, pad the input region sequence up and downstream
    return Coordinate(chrom, start - padding, stop + padding, padding)
