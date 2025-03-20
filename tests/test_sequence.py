from hapsolver.sequence import Sequence, SequenceIterator, Fasta
from hapsolver.coordinate import Coordinate

import pytest
import pysam

import os

# ====== Sequence class test

TEST_CASES = [
    ("ACGT", "ACGT", 4),  # ID: valid_sequence
    ("ACGTN", "ACGTN", 5),  # ID: valid_sequence_with_n
    ("A", "A", 1),  # ID: single_nucleotide
    ("", "", 0),  # ID: empty_sequence
    ("AaCgTn", "AACGTN", 6),  # ID: mixed_case
]

ERROR_CASES = [
    ("ACGU", ValueError),  # ID: invalid_character_u
    ("ACGT!", ValueError),  # ID: invalid_character_!
    ("123", ValueError),  # ID: invalid_numeric
]


@pytest.mark.parametrize(
    "sequence, expected_sequence, expected_len",
    TEST_CASES,
    ids=[case[2] for case in TEST_CASES],
)
def test_valid_sequence(sequence, expected_sequence, expected_len):
    # Act
    seq = Sequence(sequence)

    # Assert
    assert str(seq) == expected_sequence
    assert len(seq) == expected_len
    assert seq.sequence == expected_sequence
    assert list(seq) == list(expected_sequence)
    for i in range(len(seq)):
        assert seq[i] == expected_sequence[i]


@pytest.mark.parametrize(
    "sequence, expected_exception",
    ERROR_CASES,
    ids=[case[1].__name__ for case in ERROR_CASES],
)
def test_invalid_sequence(sequence, expected_exception):
    # Act & Assert
    with pytest.raises(expected_exception):
        Sequence(sequence)


@pytest.mark.parametrize(
    "sequence, index, expected_nucleotide",
    [
        ("ACGT", 0, "A"),  # ID: first_nucleotide
        ("ACGT", 3, "T"),  # ID: last_nucleotide
        ("ACGTN", 4, "N"),  # ID: n_nucleotide
        ("A", 0, "A"),  # ID: single_nucleotide_access
    ],
    ids=[
        "first_nucleotide",
        "last_nucleotide",
        "n_nucleotide",
        "single_nucleotide_access",
    ],
)
def test_getitem(sequence, index, expected_nucleotide):
    # Arrange
    seq = Sequence(sequence)

    # Act & Assert
    assert seq[index] == expected_nucleotide


@pytest.mark.parametrize(
    "sequence, slice_start, slice_stop, expected_slice",
    [
        ("ACGT", 0, 2, ["A", "C"]),  # ID: slice_start_to_middle
        ("ACGT", 2, 4, ["G", "T"]),  # ID: slice_middle_to_end
        ("ACGTN", 1, 4, ["C", "G", "T"]),  # ID: slice_with_n
        ("ACGT", 0, 4, ["A", "C", "G", "T"]),  # ID: slice_whole_sequence
    ],
    ids=[
        "slice_start_to_middle",
        "slice_middle_to_end",
        "slice_with_n",
        "slice_whole_sequence",
    ],
)
def test_getitem_slice(sequence, slice_start, slice_stop, expected_slice):
    # Arrange
    seq = Sequence(sequence)

    # Act & Assert
    assert seq[slice_start:slice_stop] == expected_slice


@pytest.mark.parametrize(
    "sequence, index",
    [
        ("ACGT", 4),  # ID: index_out_of_range_positive
        ("ACGT", -5),  # ID: index_out_of_range_negative
        ("A", 1),  # ID: index_out_of_range_single_nucleotide
        ("", 0),  # ID: index_out_of_range_empty_sequence
    ],
    ids=[
        "index_out_of_range_positive",
        "index_out_of_range_negative",
        "index_out_of_range_single_nucleotide",
        "index_out_of_range_empty_sequence",
    ],
)
def test_getitem_index_error(sequence, index):
    # Arrange
    seq = Sequence(sequence)

    # Act & Assert
    with pytest.raises(IndexError):
        seq[index]


# ====== SequenceIterator class test

TEST_CASES = [
    ("ACGT", ["A", "C", "G", "T"]),  # ID: valid_sequence
    ("A", ["A"]),  # ID: single_nucleotide
    ("", []),  # ID: empty_sequence
    ("ACGTN", ["A", "C", "G", "T", "N"]),  # ID: sequence_with_n
]

ERROR_CASES = [
    (object(), AttributeError),  # ID: invalid_sequence_object
    ({"invalid": "sequence"}, AttributeError),  # ID: invalid_sequence_dict
]


@pytest.mark.parametrize(
    "sequence_string, expected_nucleotides",
    TEST_CASES,
    ids=[str(case[1]) for case in TEST_CASES],
)
def test_sequence_iterator(sequence_string, expected_nucleotides):

    # Arrange
    sequence = Sequence(sequence_string)
    iterator = SequenceIterator(sequence)

    # Act & Assert
    for expected_nucleotide in expected_nucleotides:
        assert next(iterator) == expected_nucleotide

    with pytest.raises(StopIteration):
        next(iterator)


@pytest.mark.parametrize(
    "invalid_sequence, expected_exception",
    ERROR_CASES,
    ids=[case[1].__name__ for case in ERROR_CASES],
)
def test_sequence_iterator_invalid_sequence(invalid_sequence, expected_exception):

    # Act & Assert
    with pytest.raises(expected_exception):
        SequenceIterator(invalid_sequence)


# ====== Fasta class test
FASTA_CONTENT = """>chr1
ACGTACGTACGT
>chr2
AGCTAGCTAGCT
"""

FAI = "fai"

TEST_CASES = [
    ("test.fasta", None, "chr1", 1, 5, "CGTA"),  # ID: fetch_valid_coordinates
    ("test.fasta", None, "chr2", 0, 4, "AGCT"),  # ID: fetch_start_at_zero
    ("test.fasta", None, "chr1", 5, 10, "CGTAC"),  # ID: fetch_middle_sequence
    ("test.fasta", None, "chr2", 8, 12, "AGCT"),  # ID: fetch_end_of_sequence
]

ERROR_CASES = [
    ("test.fasta", None, "chr3", 1, 5, ValueError),  # ID: fetch_invalid_contig
    ("test.fasta", None, "chr1", 15, 20, ValueError),  # ID: fetch_out_of_range
    ("test.fasta", None, "chr2", -1, 5, ValueError),  # ID: fetch_negative_start
]

INDEX_TEST_CASES = [
    ("test.fasta", None),  # ID: index_no_existing_index
    ("test.fasta", "test.fasta.fai"),  # ID: index_existing_index
]


@pytest.fixture(scope="function")
def fasta_file(tmp_path):
    fasta_path = tmp_path / "test.fasta"
    with open(fasta_path, "w") as f:
        f.write(FASTA_CONTENT)
    return str(fasta_path)


@pytest.fixture(scope="function")
def fasta_index(fasta_file, tmp_path):
    index_path = tmp_path / "test.fasta.fai"
    pysam.faidx(str(fasta_file))
    return str(index_path)


@pytest.mark.parametrize(
    "fasta_file, faidx, contig, start, stop, expected_sequence",
    TEST_CASES,
    indirect=["fasta_file"],
    ids=[case[5] for case in TEST_CASES],
)
def test_fetch_valid(fasta_file, faidx, contig, start, stop, expected_sequence):
    # Arrange
    fasta = Fasta(fasta_file, faidx)
    coord = Coordinate(contig=contig, start=start, stop=stop, padding=0)

    # Act
    sequence = fasta.fetch(coord)

    # Assert
    assert str(sequence) == expected_sequence


# @pytest.mark.parametrize("fasta_file, faidx, contig, start, stop, expected_exception", ERROR_CASES, indirect=["fasta_file"], ids=[case[5].__name__ for case in ERROR_CASES])
# def test_fetch_invalid(fasta_file, faidx, contig, start, stop, expected_exception):
#     # Arrange
#     fasta = Fasta(fasta_file, faidx)
#     coord = Coordinate(contig=contig, start=start, stop=stop, padding=0)

#     # Act & Assert
#     with pytest.raises(expected_exception):
#         fasta.fetch(coord)


def test_repr(fasta_file):
    # Arrange
    fasta = Fasta(fasta_file)

    # Act
    representation = repr(fasta)

    # Assert
    assert representation == "<Fasta object>"
