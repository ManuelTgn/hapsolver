from hapsolver import Region, RegionList, Coordinate, Sequence

import pytest


# Parametrized tests for Region initialization, length, string representation, and representation
@pytest.mark.parametrize(
    "sequence, coord, expected_len, expected_str, expected_repr",
    [
        (Sequence("ATGC"), Coordinate("chr1", 10, 13, 0), 4, ">chr1:10-13\nATGC", "<Region object; region=chr1:10-13>"),  # Happy path
        (Sequence(""), Coordinate("chrX", 100, 100, 0), 0, ">chrX:100-100\n", "<Region object; region=chrX:100-100>"),  # Empty sequence
    ],
    ids=["non_empty_sequence", "empty_sequence"]
)
def test_region_init_len_str_repr(sequence: Sequence, coord: Coordinate, expected_len: int, expected_str: str, expected_repr: str) -> None:
    # Act
    region = Region(sequence, coord)

    # Assert
    assert len(region) == expected_len
    assert str(region) == expected_str
    assert repr(region) == expected_repr


# Parametrized tests for contain method
@pytest.mark.parametrize(
    "region1_coord, region2_coord, expected_contain",
    [
        (Coordinate("chr1", 10, 20, 0), Coordinate("chr1", 12, 18, 0), True),  # Happy path, full containment
        (Coordinate("chr1", 10, 20, 0), Coordinate("chr1", 10, 20, 0), True),  # Edge case, identical regions
        (Coordinate("chr1", 10, 20, 0), Coordinate("chr1", 12, 22, 0), False),  # Edge case, partial overlap, extending beyond end
        (Coordinate("chr1", 10, 20, 0), Coordinate("chr1", 8, 18, 0), False), # Edge case, partial overlap, extending before start
        (Coordinate("chr1", 10, 20, 0), Coordinate("chr2", 12, 18, 0), False),  # Edge case, different contigs
    ],
    ids=["full_containment", "identical_regions", "extends_beyond_end", "extends_before_start", "different_contigs"]
)
def test_contain(region1_coord: Coordinate, region2_coord: Coordinate, expected_contain: bool) -> None:
    # Arrange
    region1 = Region(Sequence("A" * (region1_coord.stop - region1_coord.start +1 )), region1_coord)
    region2 = Region(Sequence("A" * (region2_coord.stop - region2_coord.start + 1)), region2_coord)

    # Act
    contain = region1.contain(region2)

    # Assert
    assert contain == expected_contain


# Parametrized tests for overlap method
@pytest.mark.parametrize(
    "region1_coord, region2_coord, expected_overlap",
    [
        (Coordinate("chr1", 10, 20, 0), Coordinate("chr1", 15, 25, 0), True),  # Happy path, partial overlap
        (Coordinate("chr1", 10, 20, 0), Coordinate("chr1", 10, 20, 0), True),  # Edge case, identical regions
        (Coordinate("chr1", 10, 20, 0), Coordinate("chr1", 21, 30, 0), False),  # Edge case, no overlap
        (Coordinate("chr1", 10, 20, 0), Coordinate("chr2", 15, 25, 0), False),  # Edge case, different contigs
        (Coordinate("chr1", 10, 20, 0), Coordinate("chr1", 5, 10, 0), True), # Edge case, region2 starts before region1 but overlaps at the start
        (Coordinate("chr1", 10, 20, 0), Coordinate("chr1", 20, 25, 0), True), # Edge case, region2 starts at the end of region1 and extends beyond
    ],
    ids=["partial_overlap", "identical_regions", "no_overlap", "different_contigs", "overlap_at_start", "overlap_at_end"]
)
def test_overlap(region1_coord: Coordinate, region2_coord: Coordinate, expected_overlap: bool) -> None:
    # Arrange
    region1 = Region(Sequence("A" * (region1_coord.stop - region1_coord.start + 1)), region1_coord)
    region2 = Region(Sequence("A" * (region2_coord.stop - region2_coord.start + 1)), region2_coord)

    # Act
    overlap = region1.overlap(region2)

    # Assert
    assert overlap == expected_overlap


# Parametrized tests for contig, start, stop, sequence properties
@pytest.mark.parametrize(
    "sequence, coord, expected_contig, expected_start, expected_stop, expected_sequence",
    [
        (Sequence("ATGC"), Coordinate("chr1", 10, 13, 0), "chr1", 10, 13, Sequence("ATGC")),  # Happy path
    ],
    ids=["basic_test"]
)
def test_properties(sequence: Sequence, coord: Coordinate, expected_contig: str, expected_start: int, expected_stop: int, expected_sequence: Sequence) -> None:
    # Arrange
    region = Region(sequence, coord)

    # Act
    contig = region.contig
    start = region.start
    stop = region.stop
    sequence = region.sequence

    # Assert
    assert contig == expected_contig
    assert start == expected_start
    assert stop == expected_stop
    assert sequence == expected_sequence

# Parametrized tests for contain method error cases
@pytest.mark.parametrize(
    "region_query, expected_error",
    [
        ("not_a_region", TypeError),  # Error case, invalid region type
        (123, TypeError),  # Error case, invalid region type
    ],
    ids=["invalid_region_type_string", "invalid_region_type_int"]
)
def test_contain_errors(region_query, expected_error) -> None:
    # Arrange
    region1 = Region(Sequence("ATGC"), Coordinate("chr1", 10, 13, 0))

    with pytest.raises(expected_error):
        # Act
        region1.contain(region_query)

# Parametrized tests for overlap method error cases
@pytest.mark.parametrize(
    "region_query, expected_error",
    [
        ("not_a_region", TypeError),  # Error case, invalid region type
        (123, TypeError),  # Error case, invalid region type
    ],
    ids=["invalid_region_type_string", "invalid_region_type_int"]
)
def test_overlap_errors(region_query, expected_error) -> None:
    # Arrange
    region1 = Region(Sequence("ATGC"), Coordinate("chr1", 10, 13, 0))

    with pytest.raises(expected_error):
        # Act
        region1.overlap(region_query)

