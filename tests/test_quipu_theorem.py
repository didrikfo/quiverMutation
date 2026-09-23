"""Inverting the quipu theorem of arXiv:2305.06642.

Theorem `thm:QuipuToAn` sends a quipu to the LNA with almost separate relations
it is derived equivalent to.  quipuForms.quipuForAlmostSeparateLNA goes the
other way, which is what lets a class be identified without any mutation search.
These check the round trip, and check the answer against the paper's own table
and against what a mutation search actually reaches.
"""

import pytest

from quivermutation import quipuForms as qf
import quivermutation as qm
from helpers import line_algebra, quiet, relation_string
from paper_classification import PAPER_CLASSES, rel_lengths


def lna_from_quipu(k, m):
    """The LNA the theorem assigns to P^(m)_(k), as (length, relLengths).

    n_i = k_0 + sum_{j=1..i} (m_{j-1} + k_j + 1), relations of length m_i + 2
    starting at k_0 and at each n_i for 1 <= i <= r.
    """
    n = [k[0]]
    for i in range(1, len(k)):
        n.append(n[-1] + m[i - 1] + k[i] + 1)
    length = n[-1]
    relLengths = [0] * (length - 2)
    for index, cord in enumerate(m):
        relLengths[n[index] - 1] = cord + 2
    return length, relLengths


PAPER_QUIPUS = [
    ((1, 2, 0, 1), (2, 1, 3)),   # the worked example, -> A_{13,(1,6,8)}^{(4,3,5)}
    ((1, 1), (2,)),              # D_5
    ((1, 1), (3,)),              # D_6
    ((1, 2), (2,)),              # E_6
    ((1, 0, 1), (1, 1)),         # D~_5
    ((1, 0, 2), (1, 1)),         # the n = 7 class the paper writes in P notation
    ((1, 1, 2), (1, 1)),
    ((1, 0, 3), (1, 1)),
    ((1, 0, 0, 1), (1, 1, 1)),
    ((2, 0, 2), (1, 1)),
    ((1, 0, 2), (1, 2)),
    ((2, 3), (2,)),
]


@pytest.mark.parametrize("k, m", PAPER_QUIPUS)
def test_quipu_to_lna_and_back(k, m):
    """Round trip: quipu -> its LNA -> back to the same quipu."""
    length, relLengths = lna_from_quipu(k, m)
    assert length == len(m) + sum(k) + sum(m)  # the paper's vertex count
    assert qf.quipuForAlmostSeparateLNA(length, relLengths) == qf.canonicalQuipuParameters(k, m)


def test_the_papers_worked_example():
    """P_(1,2,0,1)^(2,1,3) is derived equivalent to A_{13,(1,6,8)}^{(4,3,5)}."""
    relLengths = [0] * 11
    for start, length in [(1, 4), (6, 3), (8, 5)]:
        relLengths[start - 1] = length
    assert qf.quipuForAlmostSeparateLNA(13, relLengths) == qf.canonicalQuipuParameters(
        (1, 2, 0, 1), (2, 1, 3)
    )


@pytest.mark.parametrize("length", sorted(PAPER_CLASSES))
def test_every_member_of_a_published_class_gives_the_same_quipu(length):
    """The theorem must be constant on each published derived equivalence class."""
    for label, entries in PAPER_CLASSES[length].items():
        quipus = set()
        for starts, lengths in entries:
            quipus.add(qf.quipuForAlmostSeparateLNA(length, rel_lengths(length, starts, lengths)))
        assert len(quipus) == 1, f"n={length} class {label} gave {quipus}"
        assert None not in quipus, f"n={length} class {label} not recognised"


@pytest.mark.parametrize("length", sorted(PAPER_CLASSES))
def test_distinct_published_classes_give_distinct_quipus(length):
    """Distinct classes must get distinct quipus -- that is the whole point.

    Unlike the Coxeter polynomial, this is a complete invariant, so it has to
    separate every pair, not just happen to.
    """
    quipus = {
        label: qf.quipuForAlmostSeparateLNA(length, rel_lengths(length, *entries[0]))
        for label, entries in PAPER_CLASSES[length].items()
    }
    assert len(set(quipus.values())) == len(quipus)


def test_relations_of_length_two_are_ignored():
    """They do not change the class, so they must not change the quipu."""
    plain = qf.quipuForAlmostSeparateLNA(6, [3, 0, 0, 0])          # A_{6,(1)}^{(3)}
    decorated = qf.quipuForAlmostSeparateLNA(6, [3, 0, 0, 2])      # the same, plus a 2
    assert plain == decorated == qf.canonicalQuipuParameters((1, 1), (3,))


def test_overlapping_relations_are_not_covered():
    """Two relations overlapping in more than one arrow are outside the theorem."""
    # Relations of 3 arrows at vertices 1 and 2 share two arrows.
    assert qf.quipuForAlmostSeparateLNA(6, [3, 3, 0, 0]) is None


@pytest.mark.slow
@pytest.mark.parametrize(
    "length, rels, depth",
    [(5, "300", 5), (5, "000", 2), (6, "3000", 5), (6, "2300", 5), (6, "3030", 7),
     (7, "30000", 5)],
)
def test_the_theorem_agrees_with_what_a_search_reaches(length, rels, depth):
    """The cheap answer and the expensive one must be the same answer.

    The theorem is O(1); reaching the same hereditary quiver by mutation costs a
    depth-5-to-7 traversal.  They are independent routes to the class, so
    agreement is a real check on both.
    """
    fromTheorem = qm.hereditaryFormFromTheorem(length, relation_string(rels))
    found = quiet(qm.hereditaryFormsReachedFrom, line_algebra(length, rels), depth)
    fromSearch = qm.formatHereditaryForms(found)
    assert fromSearch, f"the search reached no hereditary quiver from A_{length}_{rels}"
    assert fromTheorem == fromSearch


@pytest.mark.slow
def test_the_theorem_answers_where_the_search_gives_up():
    """A_{7,(2,4)}^{(3,3)} reaches no relation-free quiver within depth 8.

    This is why the theorem matters: the class is perfectly well determined, but
    finding its hereditary representative by mutation is out of reach at any
    depth the search can afford.  The theorem names it immediately.
    """
    # Deduplicated, which reaches exactly the hereditary forms the plain walk
    # does (E-043, `test_fingerprint`) in a sixth of the time at this depth.
    from quivermutation import fingerprint, search

    found = []
    search.mutationSearchDepthFirst(line_algebra(7, "03030"), 8, [], 'theorem',
                                    printOutput = False, collectedHereditary = found,
                                    visited = fingerprint.Visited())
    assert found == []
    # The paper's table puts A_{7,(2,4)}^{(3,3)} in the class of P_(1,0,2)^(1,1).
    assert qm.hereditaryFormFromTheorem(7, relation_string("03030")) == qf.formatQuipu(
        qf.canonicalQuipuParameters((1, 0, 2), (1, 1))
    ) == "P^(1,2)_(1,0,1)"
