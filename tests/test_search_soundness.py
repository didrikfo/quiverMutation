"""The search must not leave the derived equivalence class it started in.

`mutationIsPossibleAtVertex` rules mutation *out*, not in, so admissibility alone
does not make a step a derived equivalence -- R-005 recorded that for rule
discovery and R-012 for the search, where it had gone unchecked. F-038 measured
what it let through, and F-039 then found that most of that measurement was of
the *key* being computed wrong rather than of the mutation being bad: a parallel
pair of arrows counted as one path, and the cheap path count did not close the
commutativity relations to a fixed point. With both corrected there is nothing
left at `n = 6` to depth 5, and the smallest surviving case is at `n = 7` and
depth 6.

These tests pin both halves: that the guard holds the invariant, and that the
surviving counterexample is still there with the guard off, so the tests keep
failing for the right reason if the guard is ever removed.
"""

import pytest

import quivermutation as qm
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


def test_the_guard_refuses_the_one_step_that_leaves_the_class():
    """The step of F-038's surviving case, asked of the guard directly.

    The guard drops a child before the search visits it, so asking a guarded
    search whether it ever visits a moved key checks the guard against itself,
    and at `n <= 7` to depth 5 there is nothing to drop anyway (F-039).  This is
    the node the slow test below walks to with the guard off: the relation dual
    of `33030` after `[4, 1, 3, 1, 3]`, still in the class, one mutation at 3
    from leaving it.  With the guard the search refuses that mutation; without
    it, it takes it.
    """
    algebra = nk.LinearNakayamaAlgebra(7, [3, 3, 0, 3, 0])
    dual = search.memberAndItsDual(7, algebra.relationString())[1]
    base = search._coxeterKeyOrNone(dual)
    walked = dual
    for vertex in [4, 1, 3, 1, 3]:
        walked = lm._quiet(qm.reducePathAlgebra,
                           lm._quiet(qm.quiverMutationAtVertex, walked, vertex))
    assert search._coxeterKeyOrNone(walked) == base

    def moved(guard):
        seen = []
        search.mutationSearchDepthFirst(
            lm._copy(walked), 1, [], 'soundness', printOutput = False,
            visitor = lambda pathAlg, path: seen.append(
                (search._coxeterKeyOrNone(pathAlg), list(path))),
            coxeterGuard = guard, baseKey = base)
        return [path for key, path in seen if key is not None and key != base]

    assert moved(guard = False) == [[3]]
    assert moved(guard = True) == []


def test_f038s_smallest_case_was_a_mis_measurement_and_now_holds():
    """`3030` by `[1, 3, 4, 1, 4]` produces a parallel pair, and that is all.

    It was the first node E-033 found with a moved key, and the move was in the
    Cartan matrix rather than in the algebra: two arrows `1 -> 6`, counted once.
    Nothing at `n = 6` to depth 5 leaves the class now, with the guard off --
    which is why the guard's own justification is the slow `n = 7` test below.
    """
    algebra = nk.LinearNakayamaAlgebra(6, [3, 0, 3, 0])
    walked = algebra
    for vertex in [1, 3, 4, 1, 4]:
        walked = lm._quiet(qm.reducePathAlgebra,
                           lm._quiet(qm.quiverMutationAtVertex, walked, vertex))
    assert walked.hasParallelArrows()
    assert search._coxeterKeyOrNone(walked) == search._coxeterKeyOrNone(algebra)

    base, seen = _keysReached(6, (3, 0, 3, 0), 5, guard = False)
    moved = [path for key, path in seen if key is not None and key != base]
    assert moved == []


@pytest.mark.slow
def test_without_the_guard_a_clean_quiver_leaves_the_class_at_n_7():
    """Acyclic, no parallel arrows, and outside the class -- the dangerous kind.

    The one case of F-038 that survives the corrected count of F-039: the cheap
    and the exact Cartan matrices agree here, so it is the algebra that changed.
    This is what the guard is for.
    """
    base, seen = _keysReached(7, (3, 3, 0, 3, 0), 6, guard = False)
    moved = [path for key, path in seen if key is not None and key != base]
    assert [4, 1, 3, 1, 3, 3] in moved



# ---------------------------------------------------------------------------
# And costs no answers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("length, depth", [
    (5, 4), pytest.param(6, 4, marks = pytest.mark.slow)])
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
