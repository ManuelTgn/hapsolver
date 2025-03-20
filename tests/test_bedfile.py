from hapsolver.bedfile import Bed
from hapsolver.coordinate import Coordinate
from hapsolver.sequence import Fasta
from hapsolver.region import RegionList, Region

import pytest

# Example BED and FASTA data for testing
BED_CONTENT = """chr1\t100\t200
chr2\t300\t400
"""
FASTA_CONTENT = """>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chr2
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
"""

TEST_CASES = [
    ("test.bed", 0, 2),  # ID: valid_bed_no_padding
    ("test.bed", 10, 2),  # ID: valid_bed_with_padding
]

ERROR_CASES_READ = [
    ("invalid_permissions.bed", 0, PermissionError),  # ID: invalid_permissions
]

ERROR_CASES_EXTRACT = [
    ("test.bed", AttributeError),  # ID: missing_coordinates
]


@pytest.fixture(scope="function")
def bed_file(tmp_path):
    bed_path = tmp_path / "test.bed"
    bed_path.write_text(BED_CONTENT)
    return str(bed_path)


@pytest.fixture(scope="function")
def fasta_file(tmp_path):
    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(FASTA_CONTENT)
    return str(fasta_path)


@pytest.mark.parametrize(
    "bed_file, padding, expected_len",
    TEST_CASES,
    indirect=["bed_file"],
    ids=[case[2] for case in TEST_CASES],
)
def test_valid_bed(bed_file, padding, expected_len):
    # Act
    bed = Bed(bed_file, padding)

    # Assert
    assert len(bed) == expected_len
    assert repr(bed) == f"<Bed object; stored regions={expected_len}>"
    for coord in bed:
        assert isinstance(coord, Coordinate)


@pytest.mark.parametrize(
    "bed_file, padding, expected_exception",
    ERROR_CASES_READ,
    indirect=["bed_file"],
    ids=[case[2].__name__ for case in ERROR_CASES_READ],
)
def test_invalid_bed_read(bed_file, padding, expected_exception, monkeypatch):
    # sourcery skip: simplify-generator
    if expected_exception == PermissionError:
        # Mock file permissions to raise PermissionError
        monkeypatch.setattr(
            "builtins.open",
            lambda *args, **kwargs: (_ for _ in []).throw(PermissionError),
        )

    # Act & Assert
    with pytest.raises(expected_exception):
        Bed(bed_file, padding)


@pytest.mark.parametrize(
    "bed_file, expected_exception",
    ERROR_CASES_EXTRACT,
    indirect=["bed_file"],
    ids=[case[1].__name__ for case in ERROR_CASES_EXTRACT],
)
def test_extract_regions_error(bed_file, expected_exception, fasta_file):
    # Arrange
    bed = Bed(bed_file, 0)
    bed._coordinates = None  # Simulate missing coordinates

    # Act & Assert
    with pytest.raises(expected_exception):
        bed.extract_regions(Fasta(fasta_file))


def test_extract_regions(bed_file, fasta_file):
    # Arrange
    bed = Bed(bed_file, 0)
    fasta = Fasta(fasta_file)

    # Act
    regions = bed.extract_regions(fasta)

    # Assert
    assert isinstance(regions, RegionList)
    assert len(regions) == len(bed)
    for region in regions:
        assert isinstance(region, Region)


@pytest.mark.parametrize(
    "bed_file, padding, index, expected_contig, expected_start, expected_stop",
    [
        ("test.bed", 0, 0, "chr1", 100, 200),  # ID: getitem_first_element
        ("test.bed", 0, 1, "chr2", 300, 400),  # ID: getitem_second_element
        ("test.bed", 10, 0, "chr1", 90, 210),  # ID: getitem_with_padding
    ],
    indirect=["bed_file"],
    ids=["getitem_first_element", "getitem_second_element", "getitem_with_padding"],
)
def test_getitem(
    bed_file, padding, index, expected_contig, expected_start, expected_stop
):
    # Arrange
    bed = Bed(bed_file, padding)

    # Act
    coordinate = bed[index]

    # Assert
    assert coordinate.contig == expected_contig
    assert coordinate.start == expected_start
    assert coordinate.stop == expected_stop
