""" """

from typing import List, Any

# define static variables shared across modules
DNA = ["A", "C", "G", "T", "N"]  # dna alphabet


# define utils functions
def flatten_list(lst: List[List[Any]]) -> List[Any]:
    """Flattens a list of lists into a single list.

    Args:
        lst: The list of lists to flatten.
    Returns:
        A new list containing all the elements of the sublists in a single flattened list.
    """
    return [e for sublist in lst for e in sublist]


def adjust_position(position: int, pivot_pos: int) -> int:
    """Adjusts a position relative to a pivot position.

    Args:
        position: The position to adjust.
        pivot_pos: The pivot position.
    Returns:
        The adjusted position.
    """
    return position - pivot_pos

def compute_indel_length(ref: str, alt: str) -> int:
    """Computes the length of an indel.

    Args:
        ref: The reference sequence.
        alt: The alternate sequence.
    Returns:
        The length of the indel.
    """
    return abs(len(ref) - len(alt))
