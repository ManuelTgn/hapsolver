from hapsolver import Node, HaplotypeGraph

from typing import List, Set, Tuple

import pytest

class VariantRecord:  # test class
    def __init__(self, id, pos, alleles, samples):
        self.id = id
        self.pos = pos
        self.alleles = alleles
        self.samples = samples


@pytest.mark.parametrize(
    "variants, samples, expected_repr",
    [
        (
            {VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[(set(), set()), ({"sample1"}, {"sample2"})])},
            {"sample1"},
            "<Node object; variants=rs1, samples={'sample1'}>",
        ),
        (
            {VariantRecord(id="rs2", pos=20, alleles={"G", "T"}, samples=[({"sample3"}, {"sample4"})])},
            {"sample3"},
            "<Node object; variants=rs2, samples={'sample3'}>",
        ),
    ],
    ids=["single_variant", "single_variant_different_alleles"],
)
def test_node_repr(variants: Set[VariantRecord], samples: Set[str], expected_repr: str):
    # Act
    node = Node(variants, samples)

    # Assert
    assert repr(node) == expected_repr


@pytest.mark.parametrize(
    "variants1, samples1, variants2, samples2, expected_common_samples",
    [
        (
            {VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[(set(), set()), ({"sample1"}, {"sample2"})])},
            {"sample1"},
            {VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[(set(), set()), ({"sample1"}, {"sample2"})])},
            {"sample1", "sample2"},
            {"sample1"},
        ),
        (
            {VariantRecord(id="rs2", pos=20, alleles={"G", "T"}, samples=[({"sample3"}, {"sample4"})])},
            {"sample3"},
            {VariantRecord(id="rs3", pos=30, alleles={"A", "C"}, samples=[({"sample5"}, set())])},
            {"sample5"},
            set(),
        ),
    ],
    ids=["common_samples", "no_common_samples"],
)
def test_node_common_samples(
    variants1: Set[VariantRecord],
    samples1: Set[str],
    variants2: Set[VariantRecord],
    samples2: Set[str],
    expected_common_samples: Set[str],
):
    # Arrange
    node1 = Node(variants1, samples1)
    node2 = Node(variants2, samples2)

    # Act
    common_samples = node1.common_samples(node2)

    # Assert
    assert common_samples == expected_common_samples


@pytest.mark.parametrize(
    "variants, samples, n, expected_counter",
    [
        (
            {VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[(set(), set()), ({"sample1"}, {"sample2"})])},
            {"sample1"},
            1,
            0,
        ),
        (
            {VariantRecord(id="rs2", pos=20, alleles={"G", "T"}, samples=[({"sample3"}, {"sample4"})])},
            {"sample3", "sample4", "sample5"},
            2,
            1,
        ),
    ],
    ids=["decrease_to_zero", "decrease_partially"],
)
def test_node_decrease_counter(variants: Set[VariantRecord], samples: Set[str], n: int, expected_counter: int):
    # Arrange
    node = Node(variants, samples)

    # Act
    node.decrease_counter(n)

    # Assert
    assert node.counter == expected_counter


@pytest.mark.parametrize(
    "variants, samples, linked_node",
    [
        (
            {VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[(set(), set()), ({"sample1"}, {"sample2"})])},
            {"sample1"},
            Node({VariantRecord(id="rs2", pos=20, alleles={"G", "T"}, samples=[({"sample3"}, {"sample4"})])}, {"sample3"}),
        )
    ],
    ids=["valid_link"],
)
def test_node_link(variants: Set[VariantRecord], samples: Set[str], linked_node: Node):
    # Arrange
    node = Node(variants, samples)

    # Act
    node.link(linked_node)

    # Assert
    assert linked_node in node.edges

@pytest.mark.parametrize(
    "variants, samples, linked_node",
    [
        (
            {VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[(set(), set()), ({"sample1"}, {"sample2"})])},
            {"sample1"},
            "invalid_node",
        )
    ],
    ids=["invalid_link"],
)
def test_node_link_invalid_type(variants: Set[VariantRecord], samples: Set[str], linked_node: str):
    # Arrange
    node = Node(variants, samples)

    # Act and Assert
    with pytest.raises(TypeError):
        node.link(linked_node)


@pytest.mark.parametrize(
    "variants, chromcopy, expected_repr",
    [
        (
            [
                VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[({"sample1"}, {"sample2"})]),
                VariantRecord(id="rs2", pos=20, alleles={"G", "T"}, samples=[({"sample1"}, {"sample2"})]),
            ],
            0,
            "<HaplotypeGraph object; nodes=2, layers=2>",
        ),
        (
            [
                VariantRecord(id="rs3", pos=30, alleles={"A", "C"}, samples=[({"sample3"}, {"sample4"})]),
            ],
            1,
            "<HaplotypeGraph object; nodes=1, layers=1>",
        ),
    ],
    ids=["multiple_variants", "single_variant"],
)
def test_haplotype_graph_repr(variants: List[VariantRecord], chromcopy: int, expected_repr: str):
    # Act
    graph = HaplotypeGraph(variants, chromcopy)

    # Assert
    assert repr(graph) == expected_repr


@pytest.mark.parametrize(
    "variants, chromcopy",
    [
        (
            [
                VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[({"sample1"}, {"sample2"})]),
            ],
            2,
        )
    ],
    ids=["invalid_chromcopy"],
)
def test_haplotype_graph_invalid_chromcopy(variants: List[VariantRecord], chromcopy: int):

    # Act and Assert
    with pytest.raises(ValueError):
        HaplotypeGraph(variants, chromcopy)


@pytest.mark.parametrize(
    "variants, samples, layer, expected_nodes, expected_layers",
    [
        (
            {VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[({"sample1"}, {"sample2"})])},
            {"sample1"},
            0,
            2,
            {0: [Node({VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[({"sample1"}, {"sample2"})])}, {"sample1"})], 1: [Node({VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[({"sample1"}, {"sample2"})])}, {"sample1"})]},
        )
    ],
    ids=["add_node"],
)
def test_haplotype_graph_add_node(variants, samples, layer, expected_nodes, expected_layers):
    # Arrange
    graph = HaplotypeGraph([VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[({"sample1"}, {"sample2"})])])

    # Act
    graph.add_node(variants, samples, layer)

    # Assert
    assert len(graph._nodes) == expected_nodes
    assert graph._layers == expected_layers



@pytest.mark.parametrize(
    "variants, chromcopy, expected_haplotypes",
    [
        (
            [
                VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[({"sample1"}, {"sample2"})]),
                VariantRecord(id="rs2", pos=20, alleles={"G", "T"}, samples=[({"sample1"}, {"sample2"})]),
            ],
            0,
            [({VariantRecord(id="rs1", pos=10, alleles={"A", "C"}, samples=[({"sample1"}, {"sample2"})]), VariantRecord(id="rs2", pos=20, alleles={"G", "T"}, samples=[({"sample1"}, {"sample2"})])}, {"sample1"})],
        ),
        (
            [
                VariantRecord(id="rs3", pos=30, alleles={"A", "C"}, samples=[({"sample3"}, {"sample4"})]),
                VariantRecord(id="rs4", pos=40, alleles={"G", "T"}, samples=[({"sample3"}, set())]),

            ],
            1,
            [({VariantRecord(id="rs3", pos=30, alleles={"A", "C"}, samples=[({"sample3"}, {"sample4"})])}, {"sample4"})],
        ),
    ],
    ids=["multiple_variants", "single_variant_different_samples"],
)
def test_haplotype_graph_retrieve_haplotypes(variants: List[VariantRecord], chromcopy: int, expected_haplotypes: List[Tuple[Set[VariantRecord], Set[str]]]):
    # Arrange
    graph = HaplotypeGraph(variants, chromcopy)

    # Act
    graph.compute_graph()
    haplotypes = graph.retrieve_haplotypes()

    # Assert
    assert haplotypes == expected_haplotypes



