"""The search must not leave the derived equivalence class it started in.

`mutationIsPossibleAtVertex` rules mutation *out*, not in, so admissibility alone
does not make a step a derived equivalence -- R-005 recorded that for rule
discovery and R-012 for the search, where it had gone unchecked. F-038 measured
what it let through: at `n = 7` and depth 6, 97 quivers in the search tree that
are acyclic, have no parallel arrows, and carry a different Coxeter polynomial,
with the search descending from every one.

These tests pin both halves: that the guard holds the invariant, and that turning
it off still reproduces the counterexamples, so the tests keep failing for the
right reason if the guard is ever removed.
"""

import pytest

from quivermutation import lnaMoves as lm
from quivermutation import nakayama as nk
from quivermutation import pathAlgebra
from quivermutation import search


def _keysReached(length, relLengths, depth, guard):
    """Every (key, path) the search visits, and the start's key."""
    algebra = nk.LinearNakayamaAlgebra(length, list(relLengths))
    seen = []
    for startAlg in search.memberAndItsDual(length, algebra.relationString()):
        base = search._coxeterKeyOrNone(startAlg)

        def visitor(pathAlg, path):
            seen.append((search._coxeterKeyOrNone(pathAlg), list(path)))

        search.mutationSearchDepthFirst(lm._copy(startAlg), depth, [], 'soundness',
                                        printOutput = False, visitor = visitor,
                                        coxeterGuard = guard)
    return base, seen


def _linesReached(length, relLengths, depth, guard):
    algebra = nk.LinearNakayamaAlgebra(length, list(relLengths))
    reached = set()
    for startAlg in search.memberAndItsDual(length, algebra.relationString()):
        collected = []
        search.mutationSearchDepthFirst(lm._copy(startAlg), depth, [], 'soundness',
                                        printOutput = False, collected = collected,
                                        coxeterGuard = guard)
        for pathAlg, _path, _numbering in collected:
            row = lm.asRelLengths(lm._copy(pathAlg), length)
            if row is not None:
                reached.add(tuple(row))
    return reached


# ---------------------------------------------------------------------------
# The guard holds the invariant
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("relLengths", [(3, 0, 3, 0), (0, 3, 3, 0), (2, 0, 0, 2)])
def test_every_quiver_the_guarded_search_reaches_keeps_the_coxeter_polynomial(relLengths):
    base, seen = _keysReached(6, relLengths, 5, guard = True)
    moved = [(key, path) for key, path in seen if key is not None and key != base]
    assert moved == []


def test_without_the_guard_the_search_leaves_the_class_at_n_6():
    """F-038's smallest case, kept so the guard is testing something."""
    base, seen = _keysReached(6, (3, 0, 3, 0), 5, guard = False)
    moved = [path for key, path in seen if key is not None and key != base]
    assert [1, 3, 4, 1, 4] in moved


@pytest.mark.slow
def test_without_the_guard_a_clean_quiver_leaves_the_class_at_n_7():
    """Acyclic, no parallel arrows, and outside the class -- the dangerous kind."""
    base, seen = _keysReached(7, (3, 3, 0, 3, 0), 6, guard = False)
    moved = [path for key, path in seen if key is not None and key != base]
    assert [4, 1, 3, 1, 3, 3] in moved


@pytest.mark.slow
@pytest.mark.parametrize("length, depth", [(6, 5), (7, 5)])
def test_the_guarded_search_never_leaves_the_class_anywhere(length, depth):
    for relLengths in nk.allRelationLengths(length):
        base, seen = _keysReached(length, relLengths, depth, guard = True)
        assert not [p for key, p in seen if key is not None and key != base], relLengths


# ---------------------------------------------------------------------------
# And costs no answers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("length, depth", [(5, 4), (6, 4)])
def test_the_guard_loses_no_line_and_gains_none(length, depth):
    """F-038: over every LNA at n = 6 and 7 to depth 5, the two agree exactly."""
    for relLengths in nk.allRelationLengths(length):
        assert (_linesReached(length, relLengths, depth, guard = True)
                == _linesReached(length, relLengths, depth, guard = False)), relLengths


def test_a_cyclic_start_disables_the_guard_rather_than_raising():
    """A cycle has no unimodular Cartan matrix, so there is no invariant to use.

    `coxeterKey` raises there, and the guard computing it at the top of the
    search turned `test_cycles.py`'s "end the branch instead of crashing" into a
    crash.  `_coxeterKeyOrNone` is what keeps that test honest.
    """
    cyclic = pathAlgebra.PathAlgebra()
    cyclic.add_arrows_from([[1, 2], [2, 3], [3, 1]])
    cyclic.add_rels_from([[[1, 2, 3]]])

    assert search._coxeterKeyOrNone(cyclic) is None
    collected = []
    search.mutationSearchDepthFirst(cyclic, 3, [], 'cyclic', printOutput = False,
                                    collected = collected)
    assert collected == []
