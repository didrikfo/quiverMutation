"""Conditional deeper probing: more depth for the branches that earn it.

A depth-bounded search gives every branch the same budget.  The quivers worth a
longer look are rare -- the ones with parallel arrows are the case this was
built for, since the procedure only learned to mutate into them in F-039 -- and
raising the depth everywhere to reach past them costs about fivefold a level.
`search.DeeperWhen` raises it where a condition holds and nowhere else.

What these pin:

* the budget, which is the only reason the deeper walk terminates: a branch
  never runs longer than `depth + budget`, and with a condition that always
  holds it runs exactly that long;
* that a grant of 0 changes nothing about what the search reaches, so the same
  object is a recorder;
* that recording and re-searching *contains* deepening in one pass, rather than
  being another way of spelling it -- every recorded firing starts a fresh
  budget where one pass spends one budget down the whole branch;
* that a condition on oriented cycles can never buy depth, because the search
  does not descend from a cyclic quiver at all;
* the spec strings a command line passes in.
"""

import networkx as nx
import pytest

from quivermutation import nakayama as nk
from quivermutation import pathAlgebra as pa
from quivermutation import search


def start(name = "3030"):
    """An LNA the search walks into a parallel-arrow quiver from."""
    return nk.LinearNakayamaAlgebra(len(name) + 2, name)


def walkRecording(algebra, depth, deeperWhen = None):
    """Every node a search reaches, as (path, signature) pairs."""
    seen = []

    def visitor(pathAlg, path):
        seen.append((list(path), signature(pathAlg)))

    search.mutationSearchDepthFirst(algebra, depth, printOutput = False,
                                    visitor = visitor, deeperWhen = deeperWhen)
    return seen


def signature(pathAlg):
    """What tells one node from another, labels and arrow names included."""
    return (tuple(sorted(pathAlg.quiver.edges(keys = True))),
            tuple(sorted(tuple(tuple(path) for path in rel) for rel in pathAlg.rels)))


def deepest(seen):
    return max(len(path) for path, _ in seen)


# -- the budget -----------------------------------------------------------

def test_a_branch_never_runs_longer_than_the_depth_plus_the_budget():
    """The one property that makes the deeper walk finite.

    A grant that renewed at every node where the condition holds would not
    terminate: a branch inside the region would refill its depth faster than it
    spent it.  With a condition that holds *everywhere* the bound is tight, so
    this measures the bound rather than merely respecting it.
    """
    always = search.DeeperWhen(lambda pathAlg: True, extraDepth = 1, budget = 3)
    assert deepest(walkRecording(start(), 2)) == 2
    assert deepest(walkRecording(start(), 2, always)) == 5


def test_the_budget_defaults_to_one_grant():
    """`extraDepth` with no budget said means: bought once, not once a node."""
    once = search.DeeperWhen(lambda pathAlg: True, extraDepth = 2)
    assert once.budget == 2
    assert deepest(walkRecording(start(), 2, once)) == 4
    # The root itself satisfies the condition, so the whole tree is one branch
    # as far as the budget is concerned and exactly one grant is ever made.
    assert once.grants == 1


def test_a_larger_budget_buys_a_second_grant():
    """`grants` counts the tree, not the branch: the root spends one and every
    child of it spends the other, so the count is one per node that bought."""
    once = search.DeeperWhen(lambda pathAlg: True, extraDepth = 1, budget = 1)
    twice = search.DeeperWhen(lambda pathAlg: True, extraDepth = 1, budget = 2)
    assert deepest(walkRecording(start(), 2, once)) == 3
    assert deepest(walkRecording(start(), 2, twice)) == 4
    assert once.grants == 1
    assert twice.grants > once.grants


def test_the_limit_caps_the_grants_in_a_run():
    """The valve for an overnight job: record everything, buy at most so much."""
    capped = search.DeeperWhen(search.hasParallelArrows, extraDepth = 2, limit = 1)
    walkRecording(start(), 5, capped)
    assert capped.grants == 1
    assert len(capped.firings) > 1, "the condition should have fired more than once"
    assert sum(firing['granted'] for firing in capped.firings) == 2


def test_a_penalty_is_not_a_grant():
    with pytest.raises(ValueError):
        search.DeeperWhen(search.hasParallelArrows, extraDepth = -1)
    with pytest.raises(ValueError):
        search.DeeperWhen(search.hasParallelArrows, extraDepth = 1, budget = -1)


# -- what it reaches ------------------------------------------------------

def test_the_probe_fires_on_the_parallel_arrows_and_walks_past_them():
    probe = search.DeeperWhen(search.hasParallelArrows, extraDepth = 2)
    plain = walkRecording(start(), 4)
    deeper = walkRecording(start(), 4, probe)
    assert probe.firings, "no parallel-arrow quiver was reached at all"
    assert probe.grants > 0
    assert deepest(plain) == 4
    assert deepest(deeper) == 6
    # Nothing the plain walk saw is lost: the grant only ever adds.
    assert {sig for _, sig in plain} <= {sig for _, sig in deeper}


def test_a_grant_of_nothing_changes_nothing():
    """Which is what makes the same object a recorder -- see `recordOnly`."""
    recorder = search.recordOnly(search.hasParallelArrows, keepQuivers = False)
    plain = walkRecording(start(), 4)
    recorded = walkRecording(start(), 4, recorder)
    assert [sig for _, sig in plain] == [sig for _, sig in recorded]
    assert recorder.firings and recorder.grants == 0


def test_recording_and_searching_again_contains_deepening_in_one_pass():
    """The two ways of asking, and where the difference between them lives.

    One pass spends a single budget down a whole branch; a recorded firing
    re-searched afterwards starts a fresh one.  That can only matter for a
    branch that *leaves* the interesting region and comes back, because a firing
    below the one that bought the depth is already inside the subtree the grant
    paid for, at exactly the remaining depth re-searching it would give.  So the
    second pass contains the first, and here -- as at every size measured, E-036
    -- the two come out equal, which is what says the one-go run is not the
    weaker of the two in practice.
    """
    depth, extra = 4, 2
    plain = {sig for _, sig in walkRecording(start(), depth)}
    onePass = {sig for _, sig in
               walkRecording(start(), depth,
                             search.DeeperWhen(search.hasParallelArrows, extraDepth = extra))}

    recorder = search.recordOnly(search.hasParallelArrows)
    twoPass = {sig for _, sig in walkRecording(start(), depth, recorder)}
    for firing in recorder.firings:
        twoPass |= {sig for _, sig in
                    walkRecording(firing['pathAlg'], firing['depth'] + extra)}

    assert onePass <= twoPass
    assert plain < onePass, "otherwise this pins nothing"
    assert onePass == twoPass


def test_a_kept_quiver_is_one_the_search_can_be_resumed_from():
    recorder = search.recordOnly(search.hasParallelArrows)
    walkRecording(start(), 4, recorder)
    firing = recorder.firings[0]
    assert firing['pathAlg'].hasParallelArrows()
    assert firing['parallelArrows'] >= 1
    assert len(firing['path']) + firing['depth'] == 4
    assert walkRecording(firing['pathAlg'], 2)


# -- the conditions -------------------------------------------------------

def test_a_condition_on_cycles_can_record_but_can_never_deepen():
    """The trap the registry documents rather than hides.

    `mutationSearchDepthFirst` does not descend from a quiver with an oriented
    cycle, so a node where this fires has no children for the extra depth to be
    spent on.  The grant is made and buys nothing.
    """
    cyclic = pa.PathAlgebra()
    cyclic.add_vertices_from([1, 2, 3])
    cyclic.add_arrows_from([[1, 2], [2, 3], [3, 1]])
    probe = search.DeeperWhen(search.hasOrientedCycle, extraDepth = 4)
    seen = walkRecording(cyclic, 2, probe)
    assert probe.grants == 1
    assert deepest(seen) == 0


def test_counting_parallel_arrows_does_not_count_a_two_cycle():
    """Two opposite arrows are not parallel, though the undirected graph merges
    them -- which is what `describeRelationFreeQuiver` counts and this does not."""
    quiver = nx.MultiDiGraph()
    quiver.add_edge(1, 2)
    quiver.add_edge(2, 1)
    assert search.countParallelArrows(quiver) == 0
    quiver.add_edge(1, 2)
    assert search.countParallelArrows(quiver) == 1


def test_no_relations_is_the_hereditary_condition():
    assert search.hasNoRelations(pa.PathAlgebra())
    assert not search.hasNoRelations(start())


def test_the_combinators_read_as_they_say():
    both = search.allOf(search.hasParallelArrows, search.hasNoRelations)
    either = search.anyOf(search.hasParallelArrows, search.hasNoRelations)
    algebra = start()
    assert not both(algebra)
    assert not either(algebra)
    assert 'hasParallelArrows' in both.__name__


def test_at_least_counts_the_extra_arrows_not_the_pairs():
    quiver = nx.MultiDiGraph()
    quiver.add_edge(1, 2)
    quiver.add_edge(1, 2)
    algebra = pa.PathAlgebra()
    algebra.quiver = quiver
    assert search.parallelArrowsAtLeast(1)(algebra)
    assert not search.parallelArrowsAtLeast(2)(algebra)


# -- the spec a command line passes ---------------------------------------

def test_the_spec_names_a_condition_and_its_numbers():
    probe = search.deeperWhenFromSpec('parallel-arrows')
    assert probe.extraDepth == 2 and probe.budget == 2 and probe.limit is None
    probe = search.deeperWhenFromSpec('parallel-arrows:3')
    assert probe.extraDepth == 3 and probe.budget == 3
    probe = search.deeperWhenFromSpec('parallel-arrows:3:6:100')
    assert (probe.extraDepth, probe.budget, probe.limit) == (3, 6, 100)
    assert probe.name == 'parallel-arrows'


def test_a_spec_that_names_nothing_says_what_there_is():
    with pytest.raises(ValueError) as raised:
        search.deeperWhenFromSpec('parallel')
    assert 'parallel-arrows' in str(raised.value)
    with pytest.raises(ValueError):
        search.deeperWhenFromSpec('parallel-arrows:deep')
    with pytest.raises(ValueError):
        search.deeperWhenFromSpec('parallel-arrows:1:2:3:4')


def test_the_summary_is_what_a_run_prints():
    probe = search.deeperWhenFromSpec('parallel-arrows:2')
    walkRecording(start(), 4, probe)
    summary = probe.summarise()
    assert summary['condition'] == 'parallel-arrows'
    assert summary['firings'] == len(probe.firings)
    assert summary['depthGranted'] == 2 * summary['grants']
    assert summary['shallowest'] <= summary['deepest']


# -- the classification pipeline ------------------------------------------

def test_the_spec_round_trips_for_a_registered_condition():
    """What a checkpoint stores has to name the same search when it is read."""
    for spec in ['parallel-arrows:2:2', 'parallel-arrows:3:6', 'no-relations:1:4:50']:
        assert search.deeperWhenFromSpec(spec).spec() == spec
    # A shorthand normalises, so one search is one entry rather than two names.
    assert search.deeperWhenFromSpec('parallel-arrows').spec() == 'parallel-arrows:2:2'
    # A condition of one's own is named but cannot be read back; the docstring
    # says so, and this is what it means.
    own = search.DeeperWhen(lambda pathAlg: True, extraDepth = 1, name = 'mine')
    assert own.spec() == 'mine:1:1'
    with pytest.raises(ValueError):
        search.deeperWhenFromSpec(own.spec())


def test_a_resume_under_a_different_condition_redoes_rather_than_skips(tmp_path, monkeypatch):
    """A probed run is not the plain run at the same depth, so its records are not.

    The table is kept either way -- a row placed is placed, whatever placed it --
    and only the records that would let a step be skipped are dropped.
    """
    import quivermutation as qm
    from helpers import quiet

    monkeypatch.chdir(tmp_path)
    quiet(qm.classifyLength, 6, 6, 6, None, False)
    fileName = "A_6_mutation_classes.csv"

    fabricated = {'named': {'someClass': 99}, 'resolved': {'other': [6, 1]}, 'condition': ''}
    qm.writeProgress(fileName, dict(fabricated))
    quiet(qm.classifyLength, 6, 6, 6, None, False, True)
    kept = qm.readProgress(fileName)
    assert kept['named'] == fabricated['named'], "a plain resume should keep them"
    assert kept['resolved'] == fabricated['resolved']

    qm.writeProgress(fileName, dict(fabricated))
    probe = search.DeeperWhen(search.hasParallelArrows, extraDepth = 2)
    quiet(qm.classifyLength, 6, 6, 6, None, False, True, deeperWhen = probe)
    dropped = qm.readProgress(fileName)
    assert dropped['named'] == {}, "the records of a plain run were reused"
    assert dropped['resolved'] == {}
    assert dropped['condition'] == probe.spec()
    # And the rows are all still placed: the table is not what was dropped.
    table = qm.mutationClassTable.MutationClassTable.fromCSV(fileName, 6)
    assert not table.unassignedRelationStrings()
