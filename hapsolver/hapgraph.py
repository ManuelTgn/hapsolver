""" """

from .variant import VariantRecord

from itertools import combinations
from typing import Set, List, Optional, Tuple

import sys


class Node:
    """Represents a node in the haplotype graph.

    A node stores a set of variants, the samples carrying those variants,
    a counter of the samples, and edges to child nodes.
    """

    def __init__(self, variants: Set[VariantRecord], samples: Set[str]) -> None:
        """Initializes a new Node object.

        Args:
            variants: The set of variants pointed to by this node.
            samples: The set of samples carrying the pointed variants.
        """
        self._variants = variants  # set of variants pointed by current node
        self._samples = samples  # set of samples carrying the pointed variants
        self._counter = len(samples)  # count samples carrying the allele
        self._edges = set()  # pointer to children nodes

    def __repr__(self) -> str:
        """Returns a string representation of the Node.

        Args:
            None
        Returns:
            A string representation of the Node.
        """
        vids = ",".join("".join(variant.id) for variant in self._variants)
        return f"<{self.__class__.__name__} object; variants={vids}, samples={self._samples}>"

    def link(self, node: "Node") -> None:
        """Links the current node to another node.

        Adds the provided node to the set of edges of the current node.
        Raises a TypeError if the provided node is not of the same class.

        Args:
            node: The node to link to.
        Returns:
            None
        """
        if not isinstance(node, self.__class__):
            raise TypeError(
                f"Cannot link objects of type {type(node).__name__} and {self.__class__.__name__} objects"
            )
        self._edges.add(node)

    def common_samples(self, node_query: "Node") -> Set[str]:
        """Returns the set of samples common to both nodes.

        Computes the intersection of the samples in the current node and the provided
            node.

        Args:
            node_query: The node to compare against.
        Returns:
            A set of common samples.
        """
        return self._samples.intersection(node_query.samples)

    def decrease_counter(self, n: int) -> None:
        """Decreases the sample counter by a given amount.

        Subtracts the provided value from the current sample counter.

        Args:
            n: The amount to decrease the counter by.
        Returns:
            None
        """
        self._counter -= n

    @property
    def edges(self) -> Set["Node"]:
        return self._edges

    @property
    def variants(self) -> Set[VariantRecord]:
        return self._variants

    @property
    def samples(self) -> Set[str]:
        return self._samples

    @property
    def counter(self) -> int:
        return self._counter


class HaplotypeGraph:
    """Represents a haplotype graph.

    The graph is constructed from a list of variants and a chromosome copy.
    """

    def __init__(
        self, variants: List[VariantRecord], chromcopy: Optional[int] = 0
    ) -> None:
        """Initializes a new HaplotypeGraph object.

        Args:
            variants: A list of VariantRecord objects.
            chromcopy: The chromosome copy to use for haplotype reconstruction.
        Raises:
            ValueError: If chromcopy is not 0 or 1.
        """
        # sourcery skip: identity-comprehension
        if chromcopy not in [0, 1]:
            raise ValueError(
                f"Forbidden chromosome copy for haplotype reconstruction ({chromcopy})"
            )
        # create starting nodes in the haplotype graph (DAG)
        self._nodes = [
            Node({variant}, variant.samples[0][chromcopy]) for variant in variants
        ]
        self._layers_num = len(variants)  # number of maximum layers in the graph
        # initialize graph layers
        self._layers = {i: [] for i in range(self._layers_num)}
        # initialize first layer with starting nodes
        self._layers[0] = [n for n in self._nodes]
        # store chromosome copy where the variant occurs
        self._chromcopy = chromcopy

    def __repr__(self) -> str:
        """Returns a string representation of the HaplotypeGraph.

        Args:
            None
        Returns:
            A string representation of the HaplotypeGraph.
        """
        return f"<{self.__class__.__name__} object; nodes={len(self._nodes)}, layers={self._layers_num}>"

    def add_node(
        self, variants: Set[VariantRecord], samples: Set[str], layer: int
    ) -> Node:
        """Adds a new node to the haplotype graph.

        Creates a new node with the given variants and samples, adds it to the graph,
        and adds it to the corresponding layer.

        Args:
            variants: The set of variants for the new node.
            samples: The set of samples for the new node.
            layer: The layer to add the node to.
        Returns:
            The newly created node.
        """
        node = Node(variants, samples)  # create new node object
        self._nodes.append(node)  # add new node to haplotype graph
        self._layers[layer].append(node)  # add node to corresponding layer
        return node

    def compute_graph(self) -> None:
        """Computes the haplotype graph.

        Iterates through the layers of the graph, connecting nodes with common samples
        and creating new nodes for the intersections. The process stops when no more
        intersections are found or all samples have been mapped.
        """
        layers = self._layers_num  # initialize var storing non-empty layers
        for layer in range(self._layers_num - 1):
            for node1, node2 in list(combinations(self._layers[layer], r=2)):
                if node1.counter == 0 and node2.counter == 0:
                    continue  # all samples already mapped
                if common_samples := node1.common_samples(node2):
                    node_int = self.add_node(
                        node1.variants.union(node2.variants), common_samples, layer + 1
                    )
                    for n in [node1, node2]:  # update graph structure
                        n.link(node_int)  # connect parents to child node
                        # record samples visited
                        n.decrease_counter(len(common_samples))
            if not self._layers[layer + 1]:  # no samples intersection found
                layers = layer + 1
                break
        self._layers = {l: self._layers[l] for l in range(layers)}
        self._layers_num = len(self._layers)  # update number of layers in graph

    def retrieve_haplotypes(self) -> List[Tuple[Set[VariantRecord], Set[str]]]:
        """Retrieves the haplotypes from the graph.

        Traverses the graph by layer from bottom to top, collecting the variants and
        samples for each unreconstructed haplotype.

        Args:
            None
        Returns:
            A list of tuples, where each tuple contains a set of variants and a set
            of samples representing a haplotype.
        """
        haplotypes = []
        reconstructed_haps = []
        # traverse the haplotype graph by layer (go from bottom)
        for layer in list(range(self._layers_num))[::-1]:
            for node in self._layers[layer]:
                if unreconstructed_haps := node.samples.difference(
                    set(reconstructed_haps)
                ):
                    haplotypes.append((node.variants, unreconstructed_haps))
                    reconstructed_haps.extend(list(unreconstructed_haps))
        return haplotypes
    

    @property
    def chromcopy(self) -> int:
        return self._chromcopy
