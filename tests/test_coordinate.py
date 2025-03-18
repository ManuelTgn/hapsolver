from hapsolver.coordinate import Coordinate

import pytest


class TestCoordinate:
    @pytest.mark.parametrize(
        "contig, start, stop, padding, expected",
        [
            ("chr1", 100, 200, 0, "chr1:100-200"),  # happy path, no padding
            ("chrX", 5000, 5500, 100, "chrX:5000-5500"),  # happy path with padding
            ("chr2", 1000, 2000, 50, "chr2:1000-2000"),  # happy path with different padding
            ("chr1", 0, 1, 0, "chr1:0-1"),  # edge case, start at 0
            ("chrM", 16569, 16570, 0, "chrM:16569-16570"),  # edge case, stop = start + 1
        ],
        ids=[
            "no_padding",
            "with_padding",
            "different_padding",
            "start_at_zero",
            "stop_equals_start_plus_one",
        ],
    )
    def test_str_happy_path(self, contig, start, stop, padding, expected):
        # Arrange

        # Act
        coordinate = Coordinate(contig, start - padding, stop + padding, padding)
        result = str(coordinate)

        # Assert
        assert result == expected

    @pytest.mark.parametrize(
        "contig, start, stop, padding, expected",
        [
            ("chr1", 100, 200, 0, "<Coordinate object; coordinate=chr1:100-200; padding=0>"),  # happy path, no padding
            ("chrX", 5000, 5500, 100, "<Coordinate object; coordinate=chrX:5000-5500; padding=100>"),  # happy path with padding
            ("chr2", 1000, 2000, 50, "<Coordinate object; coordinate=chr2:1000-2000; padding=50>"),  # happy path with different padding
            ("chr1", 0, 1, 0, "<Coordinate object; coordinate=chr1:0-1; padding=0>"),  # edge case, start at 0
            ("chrM", 16569, 16570, 0, "<Coordinate object; coordinate=chrM:16569-16570; padding=0>"),  # edge case, stop = start + 1
        ],
        ids=[
            "no_padding",
            "with_padding",
            "different_padding",
            "start_at_zero",
            "stop_equals_start_plus_one",
        ],
    )
    def test_repr_happy_path(self, contig, start, stop, padding, expected):
        # Arrange

        # Act
        coordinate = Coordinate(contig, start - padding, stop + padding, padding)
        result = repr(coordinate)

        # Assert
        assert result == expected

    @pytest.mark.parametrize(
        "contig, start, stop, padding",
        [
            ("chr1", 100, 200, 0),  # happy path, no padding
            ("chrX", 5000, 5500, 100),  # happy path with padding
        ],
        ids=["no_padding", "with_padding"],
    )
    def test_properties(self, contig, start, stop, padding):
        # Arrange

        # Act
        coordinate = Coordinate(contig, start - padding, stop + padding, padding)

        # Assert
        assert coordinate.contig == contig
        assert coordinate.start == start - padding
        assert coordinate.stop == stop + padding


    @pytest.mark.parametrize(
        "start, stop, padding",
        [
            (200, 100, 0),  # error case, start > stop with no padding
            (5500, 5000, 100),  # error case, start > stop with padding
        ],
        ids=["start_greater_than_stop_no_padding", "start_greater_than_stop_with_padding"],
    )
    def test_init_start_greater_than_stop(self, start, stop, padding):
        # Act & Assert
        with pytest.raises(ValueError, match="Stop < start coordinate"):
            Coordinate("chr1", start - padding, stop + padding, padding)

