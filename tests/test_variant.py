from hapsolver.variant import (
    _assign_vtype,
    _compute_id,
    _adjust_multiallelic,
    _genotypes_to_samples,
)
from hapsolver import VariantRecord, VCF, Coordinate

from pysam import TabixFile

import warnings
import pytest
import os

warnings.filterwarnings("ignore")

VTYPES = ["snp", "indel"]  # variant types


@pytest.mark.parametrize(
    "ref, alt, expected",
    [
        ("A", "T", "snp"),  # Simple SNP
        ("AG", "A", "indel"),  # Deletion
        ("A", "AG", "indel"),  # Insertion
        ("ATC", "GTC", "snp"),  # MNP
        ("ATCG", "AT", "indel"),  # Larger deletion
        ("A", "ATCG", "indel"),  # Larger insertion
    ],
    ids=[
        "simple_snp",
        "deletion",
        "insertion",
        "mnp",
        "larger_deletion",
        "larger_insertion",
    ],
)
def test_assign_vtype(ref, alt, expected):
    # Act
    result = _assign_vtype(ref, alt)
    # Assert
    assert result == expected


@pytest.mark.parametrize(
    "chrom, pos, ref, alt, expected",
    [
        ("chr1", 100, "A", "T", "chr1-100-A/T"),  # Simple SNP
        ("chrX", 5000, "AG", "A", "chrX-5000-AG/A"),  # Deletion
        ("chr17", 100000, "T", "TA", "chr17-100000-T/TA"),  # Insertion
        ("chrM", 16000, "C", "G", "chrM-16000-C/G"),  # SNP on chrM
    ],
    ids=["simple_snp", "deletion", "insertion", "snp_chrM"],
)
def test_compute_id(chrom, pos, ref, alt, expected):
    # Act
    result = _compute_id(chrom, pos, ref, alt)
    # Assert
    assert result == expected


@pytest.mark.parametrize(
    "ref, alt, pos, expected_ref, expected_alt, expected_pos",
    [
        ("AAA", "AAT", 100, "A", "T", 102),  # Simple SNP
        ("AAAAG", "AAAA", 100, "AG", "A", 103),  # Deletion
        ("AAAAA", "AAAAAG", 100, "A", "AG", 104),  # Insertion
        ("TCTGATCG", "TCTGAT", 100, "TCG", "T", 105),  # Larger deletion
        ("TCTGA", "TCTGATCG", 100, "A", "ATCG", 104),  # Larger insertion
    ],
    ids=["simple_snp", "deletion", "insertion", "larger_deletion", "larger_insertion"],
)
def test_adjust_multiallelic(ref, alt, pos, expected_ref, expected_alt, expected_pos):
    # Act
    result_ref, result_alt, result_pos = _adjust_multiallelic(ref, alt, pos)
    # Assert
    assert (result_ref, result_alt, result_pos) == (
        expected_ref,
        expected_alt,
        expected_pos,
    )


@pytest.mark.parametrize(
    "genotypes, samples, allelesnum, phased, expected",
    [
        (
            ["1|0", "0|1"],
            ["sample1", "sample2"],
            1,
            True,
            [({"sample1"}, {"sample2"})],
        ),  # Phased het
        (
            ["1/0", "0/1"],
            ["sample1", "sample2"],
            1,
            False,
            [({"sample1", "sample2"}, set())],
        ),  # Unphased het
        (
            ["1|0", "0|0"],
            ["sample1", "sample2"],
            1,
            True,
            [({"sample1"}, set())],
        ),  # Phased with one hom ref
        (
            ["1/0", "0/0"],
            ["sample1", "sample2"],
            1,
            False,
            [({"sample1"}, set())],
        ),  # Unphased with one hom ref
        (
            ["0|0", "0|0"],
            ["sample1", "sample2"],
            1,
            True,
            [(set(), set())],
        ),  # Phased all hom ref
        (
            ["0/0", "0/0"],
            ["sample1", "sample2"],
            1,
            False,
            [(set(), set())],
        ),  # Unphased all hom ref
        (
            ["1|1", "1|0", "0|1"],
            ["s1", "s2", "s3"],
            1,
            True,
            [({"s1", "s2"}, {"s1", "s3"})],
        ),  # Phased all genotypes
        (
            ["1/1", "1/0", "0/1"],
            ["s1", "s2", "s3"],
            1,
            False,
            [({"s1", "s2", "s3"}, set())],
        ),  # Unphased all genotypes
        (
            ["1|0", "2|1", "0|2"],
            ["s1", "s2", "s3"],
            2,
            True,
            [({"s1"}, {"s2"}), ({"s2"}, {"s3"})],
        ),  # Phased multiallelic
        (
            ["1/0", "2/1", "0/2"],
            ["s1", "s2", "s3"],
            2,
            False,
            [({"s1", "s2"}, set()), ({"s3", "s2"}, set())],
        ),  # Unphased multiallelic
        (["./.", "./0"], ["s1", "s2"], 1, False, [(set(), set())]),  # Missing genotypes
    ],
    ids=[
        "phased_het",
        "unphased_het",
        "phased_homref",
        "unphased_homref",
        "phased_all_homref",
        "unphased_all_homref",
        "phased_all_genotypes",
        "unphased_all_genotypes",
        "phased_multiallelic",
        "unphased_multiallelic",
        "missing_genotypes",
    ],
)
def test_genotypes_to_samples(genotypes, samples, allelesnum, phased, expected):
    # Act
    result = _genotypes_to_samples(genotypes, samples, allelesnum, phased)
    # Assert
    assert result == expected


class TestVariantRecord:
    @pytest.mark.parametrize(
        "variant, samples, phased, chrom, position, ref, alt, allelesnum, vtype, filter, vid, samples_result",
        [
            (
                ["chr1", "100", ".", "A", "T", ".", "PASS", ".", ".", "1|0", "0|1"],
                ["sample1", "sample2"],
                True,
                "chr1",
                100,
                "A",
                ["T"],
                1,
                ["snp"],
                "PASS",
                ["chr1-100-A/T"],
                [({"sample1"}, {"sample2"})],
            ),  # Simple SNP
            (
                [
                    "chrX",
                    "5000",
                    ".",
                    "AG",
                    "A,AC",
                    ".",
                    "PASS",
                    ".",
                    ".",
                    "1/0",
                    "0/2",
                ],
                ["sample1", "sample2"],
                False,
                "chrX",
                5000,
                "AG",
                ["A", "AC"],
                2,
                ["indel", "snp"],
                "PASS",
                ["chrX-5000-AG/A", "chrX-5000-AG/AC"],
                [({"sample1"}, set()), ({"sample2"}, set())],
            ),  # Deletion
        ],
        ids=["simple_snp", "deletion"],
    )
    def test_read_vcf_line(
        self,
        variant,
        samples,
        phased,
        chrom,
        position,
        ref,
        alt,
        allelesnum,
        vtype,
        filter,
        vid,
        samples_result,
    ):

        # Act
        record = VariantRecord()
        record.read_vcf_line(variant, samples, phased)

        # Assert
        assert record.contig == chrom
        assert record.position == position
        assert record.ref == ref
        assert record.alt == alt
        assert record.allelesnum == allelesnum
        assert record.vtype == vtype
        assert record.filter == filter
        assert record.id == vid
        assert record.samples == samples_result

    @pytest.mark.parametrize(
        "altalleles, expected",
        [
            ("A,T,C", ["A", "T", "C"]),  # Multiple alleles
            ("A", ["A"]),  # Single allele
            ("", []),  # Empty string
        ],
        ids=["multiple_alleles", "single_allele", "empty_string"],
    )
    def test_retrieve_alt_alleles(self, altalleles, expected):
        # Act
        record = VariantRecord()
        result = record._retrieve_alt_alleles(altalleles)
        # Assert
        assert result == expected

    def test_retrieve_alt_alleles_error(self):
        # Act
        record = VariantRecord()
        # Assert
        with pytest.raises(AttributeError):
            record._retrieve_alt_alleles(123)  # Type error

    @pytest.mark.parametrize(
        "chrom, position, ref, alt, expected",
        [
            (
                "chr1",
                100,
                "A",
                ["T"],
                '<VariantRecord object; variant="chr1 100 A T">',
            ),  # Simple SNP
            (
                "chrX",
                5000,
                "AG",
                ["A"],
                '<VariantRecord object; variant="chrX 5000 AG A">',
            ),  # Deletion
        ],
        ids=["simple_snp", "deletion"],
    )
    def test_repr(self, chrom, position, ref, alt, expected):
        # Arrange
        record = VariantRecord()
        record._chrom = chrom
        record._position = position
        record._ref = ref
        record._alt = alt

        # Act
        result = repr(record)

        # Assert
        assert result == expected

    @pytest.mark.parametrize(
        "chrom, position, ref, alt, expected",
        [
            ("chr1", 100, "A", ["T"], "chr1\t100\tA\tT"),  # Simple SNP
            ("chrX", 5000, "AG", ["A"], "chrX\t5000\tAG\tA"),  # Deletion
        ],
        ids=["simple_snp", "deletion"],
    )
    def test_str(self, chrom, position, ref, alt, expected):
        # Arrange
        record = VariantRecord()
        record._chrom = chrom
        record._position = position
        record._ref = ref
        record._alt = alt

        # Act
        result = str(record)

        # Assert
        assert result == expected

    @pytest.mark.parametrize(
        "vtype, expected_altalleles",
        [
            ("snp", ["T"]),  # Select SNPs
            ("indel", ["A"]),  # Select indels
        ],
        ids=["select_snps", "select_indels"],
    )
    def test_get_altalleles(self, vtype, expected_altalleles):
        # Arrange
        record = VariantRecord()
        record._alt = ["T", "A"]
        record._vtype = ["snp", "indel"]

        # Act
        result = record.get_altalleles(vtype)

        # Assert
        assert result == expected_altalleles

    def test_eq(self):
        # Arrange
        record1 = VariantRecord()
        record1._chrom = "chr1"
        record1._position = 100
        record1._ref = "A"
        record1._alt = ["T"]

        record2 = VariantRecord()
        record2._chrom = "chr1"
        record2._position = 100
        record2._ref = "A"
        record2._alt = ["T"]

        # Act
        result = record1 == record2

        # Assert
        assert result

    def test_not_eq(self):
        # Arrange
        record1 = VariantRecord()
        record1._chrom = "chr1"
        record1._position = 100
        record1._ref = "A"
        record1._alt = ["T"]

        record2 = VariantRecord()
        record2._chrom = "chr2"
        record2._position = 100
        record2._ref = "A"
        record2._alt = ["T"]

        # Act
        result = record1 == record2

        # Assert
        assert not result

    def test_eq_error(self):
        # Arrange
        record1 = VariantRecord()
        record1._chrom = "chr1"
        record1._position = 100
        record1._ref = "A"
        record1._alt = ["T"]

        record2 = VariantRecord()

        # Assert
        with pytest.raises(AttributeError):
            record1 == record2

    @pytest.mark.parametrize(
        "i, expected_pos, expected_ref, expected_alt",
        [(0, 100, "A", ["T"]), (1, 100, "A", ["G"])],
        ids=["first_allele", "second_allele"],
    )
    def test_copy(self, i, expected_pos, expected_ref, expected_alt):
        # Arrange
        record = VariantRecord()
        record._chrom = "chr1"
        record._position = 100
        record._ref = "A"
        record._alt = ["T", "G"]
        record._allelesnum = 2
        record._vtype = ["snp", "snp"]
        record._filter = "PASS"
        record._vid = ["chr1-100-A/T", "chr1-101-A/G"]
        record._samples = [({"sample1"}, {"sample2"}), ({"sample3"}, {"sample4"})]

        # Act
        copy = record._copy(i)

        # Assert
        assert copy.contig == "chr1"
        assert copy.position == expected_pos
        assert copy.ref == expected_ref
        assert copy.alt == expected_alt
        assert copy.allelesnum == 1
        assert copy.filter == "PASS"

    @pytest.mark.parametrize(
        "vtype, expected_records",
        [
            ("snp", [VariantRecord(), VariantRecord()]),  # Two SNPs
            ("indel", []),  # No indels
        ],
        ids=["two_snps", "no_indels"],
    )
    def test_split(self, vtype, expected_records):
        # Arrange
        record = VariantRecord()
        record._vtype = ["snp", "snp"]
        record._ref = "A"
        record._alt = "C,G"
        record._position = 100
        record._chrom = "chr2"
        record._filter = "PASS"
        record._vid = ["chr2-100-A/C", "chr2-100-A/G"]
        record._samples = [({"s1"}, {}), ({}, {"s2"})]
        # Act
        split_records = record.split(vtype)

        # Assert
        assert len(split_records) == len(expected_records)


# Arrange
TEST_VCF = "test-data/test.vcf.gz"
TEST_VCF_INDEX = "test-data/test.vcf.gz.tbi"


# Fixtures for setting up and tearing down test data
@pytest.fixture(scope="module")
def vcf_obj():
    return VCF(TEST_VCF, TEST_VCF_INDEX)


@pytest.fixture(scope="module")
def coordinate_obj():
    return Coordinate("chr1", 10000, 20000)


# Parametrized tests for VCF initialization
@pytest.mark.parametrize(
    "fname, vcfidx, expected_contig, expected_samples, expected_phased",
    [
        (
            TEST_VCF,
            TEST_VCF_INDEX,
            "chrx",
            ["SAMPLE1", "SAMPLE2", "SAMPLE3"],
            True,
        ),  # Happy path
        (
            TEST_VCF,
            "",
            "chrx",
            ["SAMPLE1", "SAMPLE2", "SAMPLE3"],
            True,
        ),  # Happy path, index auto-detected
    ],
    ids=["with_index", "without_index"],
)
def test_vcf_init(fname, vcfidx, expected_contig, expected_samples, expected_phased):

    # Act
    vcf = VCF(fname, vcfidx)

    # Assert
    assert vcf.contig == expected_contig
    assert vcf._samples == expected_samples  # Accessing protected member for testing
    assert vcf.phased == expected_phased


# Parametrized tests for VCF initialization error cases
@pytest.mark.parametrize(
    "fname, vcfidx, expected_error",
    [
        ("non_existent.vcf.gz", "", FileNotFoundError),  # Non-existent VCF
        (TEST_VCF, "non_existent.tbi", FileNotFoundError),  # Non-existent index
    ],
    ids=["non_existent_vcf", "non_existent_index"],
)
def test_vcf_init_errors(fname, vcfidx, expected_error):
    with pytest.raises(expected_error):
        VCF(fname, vcfidx)


# Parametrized tests for _search_index method
@pytest.mark.parametrize(
    "vcfidx, expected_index",
    [
        (TEST_VCF_INDEX, TEST_VCF_INDEX),  # Happy path, index provided
        ("", TEST_VCF_INDEX),  # Happy path, index auto-detected
    ],
    ids=["with_index", "without_index"],
)
def test_search_index(vcf_obj, vcfidx, expected_index):

    # Act
    index = vcf_obj._search_index(vcfidx)  # Accessing protected member for testing

    # Assert
    assert index == expected_index


# Parametrized tests for _search_index method error cases
@pytest.mark.parametrize(
    "vcfidx, expected_error",
    [
        ("non_existent.tbi", FileNotFoundError),  # Non-existent index
    ],
    ids=["non_existent_index"],
)
def test_search_index_errors(vcf_obj, vcfidx, expected_error):
    with pytest.raises(expected_error):
        vcf_obj._search_index(vcfidx)  # Accessing protected member for testing


# Parametrized tests for index_vcf method
@pytest.mark.parametrize(
    "vcf_file",
    [
        (TEST_VCF),  # Happy path
    ],
    ids=["existing_vcf"],
)
def test_index_vcf(vcf_file):
    # Arrange
    vcf = VCF(vcf_file)
    os.remove(TEST_VCF_INDEX)  # Remove existing index for testing

    # Act
    vcf.index_vcf(pytest=True)

    # Assert
    assert os.path.exists(TEST_VCF_INDEX)


# Parametrized tests for fetch method
@pytest.mark.parametrize(
    "coordinate, expected_num_variants",
    [
        (Coordinate("chrx", 300, 400, 0), 5),  # Happy path, variants within range
        (
            Coordinate("x", 300, 400, 0),
            5,
        ),  # Happy path, variants within range, different chr prefix
        (Coordinate("chrx", 1, 300, 0), 0),  # Edge case, no variants within range
        (Coordinate("chrx", 576, 576, 0), 0),  # Edge case, start == stop
    ],
    ids=[
        "variants_within_range",
        "variants_within_range_different_chr_prefix",
        "no_variants",
        "start_equals_stop",
    ],
)
def test_fetch(vcf_obj, coordinate, expected_num_variants):
    # Act
    variants = vcf_obj.fetch(coordinate)

    # Assert
    assert len(variants) == expected_num_variants


# Parametrized tests for fetch method error cases
@pytest.mark.parametrize(
    "coordinate, expected_error",
    [
        (
            Coordinate("chr2", 10000, 20000, 0),
            ValueError,
        ),  # Error case, mismatched contig
    ],
    ids=["mismatched_contig"],
)
def test_fetch_errors(vcf_obj, coordinate, expected_error):
    with pytest.raises(expected_error):
        vcf_obj.fetch(coordinate)


# Parametrized tests for _is_phased method
@pytest.mark.parametrize(
    "vcf_file, expected_phased",
    [
        (TEST_VCF, True),  # Happy path, phased VCF
    ],
    ids=["phased_vcf"],
)
def test_is_phased(vcf_file, expected_phased):
    # Arrange
    vcf = VCF(vcf_file)

    # Act
    vcf._is_phased()  # Accessing protected member for testing

    # Assert
    assert vcf.phased == expected_phased


# Parametrized tests for contig property
@pytest.mark.parametrize(
    "vcf_file, expected_contig",
    [
        (TEST_VCF, "chrx"),  # Happy path
    ],
    ids=["chr1"],
)
def test_contig_property(vcf_file, expected_contig):
    # Arrange
    vcf = VCF(vcf_file)

    # Act
    contig = vcf.contig

    # Assert
    assert contig == expected_contig


# Parametrized tests for phased property
@pytest.mark.parametrize(
    "vcf_file, expected_phased",
    [
        (TEST_VCF, True),  # Happy path, phased VCF
    ],
    ids=["phased_vcf"],
)
def test_phased_property(vcf_file, expected_phased):
    # Arrange
    vcf = VCF(vcf_file)

    # Act
    phased = vcf.phased

    # Assert
    assert phased == expected_phased
