"""`merges.searchFrom(..., witness)` is opt-in and its paths replay to the rows they name."""
import copy
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import merges
from quivermutation import lnaMoves as lm
from quivermutation import mutation, nakayama as nk, reduction, search

LENGTH, MEMBER, DEPTH = 5, (0, 0, 0), 3


def replay(row, witness):
    """Apply the stored mutation vertices to the member (or its dual) and return the row reached."""
    relationString = nk.LinearNakayamaAlgebra(LENGTH, list(row)).relationString()
    alg = search.memberAndItsDual(LENGTH, relationString)[witness['start']]
    for v in witness['path']:
        assert mutation.mutationIsPossibleAtVertex(alg, v)
        alg = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(copy.deepcopy(alg), v))
    return alg


def test_default_has_no_witnesses_and_witness_does_not_change_the_search():
    off = merges.searchFrom((LENGTH, MEMBER, DEPTH, '', True))
    on = merges.searchFrom((LENGTH, MEMBER, DEPTH, '', True, True))
    assert off[-1] is None
    assert off[2] == on[2] and off[5] == on[5]        # same reached rows, same walk counts
    assert on[-1] is not None and set(on[-1]) == set(on[2])


def test_witness_paths_replay_to_their_rows():
    on = merges.searchFrom((LENGTH, MEMBER, DEPTH, '', True, True))
    nonTrivial = 0
    for row, w in on[-1].items():
        assert len(w['path']) <= DEPTH
        assert lm.asRelLengths(lm._copy(replay(MEMBER, w)), LENGTH) == list(row)
        nonTrivial += bool(w['path'])
    assert nonTrivial > 0
