"""This module defines the Coordinate class.

The Coordinate class represents a genomic interval. It is used when a BED file is
provided to select genomic regions for haplotype reconstruction.
"""


class Coordinate:
    """Genomic coordinate interval.

    Used when a bed file is provided in the input to select the
    genomic regions to consider when reconstructing haplotypes.

    Attributes:
        contig: Contig name.
        start: Start coordinate, including padding.
        stop: Stop coordinate, including padding.
        padding: Region padding length.
    """

    def __init__(self, contig: str, start: int, stop: int, padding: int) -> None:
        """Initialize a genomic coordinate interval.

        Used when a bed file is provided in the input to select the
        genomic regions to consider when reconstructing haplotypes.

        Attributes:
            contig: Contig name.
            start: Start coordinate, including padding.
            stop: Stop coordinate, including padding.
            padding: Region padding length.
        """
        # initialize a genomic coordinate interval
        # used when a bed file is provided in the input to select the
        # genomic regions to consider when reconstructing haplotypes
        if stop < start:
            raise ValueError("Stop < start coordinate")
        self._contig = contig  # set contig name
        self._start = start  # set start coordinate
        self._stop = stop  # set stop coordinate
        self._padding = padding  # set region padding length

    def __repr__(self) -> str:
        """Return a string representation of the Coordinate object.

        The representation includes the class name, the coordinate (contig:start-stop),
        and the padding.

        Returns:
            A string representation of the Coordinate object.
        """
        # retrieve original start and stop
        start = self._start + self._padding
        stop = self._stop - self._padding
        return f"<{self.__class__.__name__} object; coordinate={self._contig}:{start}-{stop}; padding={self._padding}>"

    def __str__(self) -> str:
        """Return a string representation of the Coordinate object.

        The representation includes the coordinate (contig:start-stop).

        Returns:
            A string representation of the Coordinate object.
        """
        # retrieve original start and stop
        start = self._start + self._padding
        stop = self._stop - self._padding
        return f"{self._contig}:{start}-{stop}"

    @property
    def contig(self) -> str:
        return self._contig

    @property
    def start(self) -> int:
        return self._start

    @property
    def stop(self) -> int:
        return self._stop
