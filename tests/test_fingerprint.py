"""Identifying an algebra the walk has met before.

The claims that matter are two.  `canonicalKey` must be **exact**: equal keys
exactly when the presentations are equal, including across the arrow-key
renaming that a quiver with parallel arrows makes ambiguous.  And a walk that
deduplicates on it must reach exactly what the plain walk reaches -- the routes
are what get deduplicated, never the destinations.
"""

import copy

import pytest

import quivermutation as qm
from quivermutation import arrowPaths
from quivermutation import fingerprint
from quivermutation import nakayama as nk
from quivermutation import procedure
from quivermutation import search

from helpers import line_algebra


# The quiver a depth-4 walk out of the n = 9 leftover orbit 3033030 reaches at
# [4, 6, 4, 7]: nine arrows, of which two are the parallel pair 7 -> 9, and a
# commutativity relation that runs along both.  It is the smallest thing in the
# repo that the old `quiverKey` had to refuse, so it is what the canonical form
# is pinned against.
PARALLEL_ARROWS = [
    (1, 2, 0), (2, 3, 0), (3, 5, 0), (4, 7, 0), (5, 8, 0),
    (6, 7, 0), (7, 9, 0), (7, 9, 1), (8, 7, 0),
]
PARALLEL_RELATIONS = [
    {((1, 2, 0), (2, 3, 0), (3, 5, 0)): 1},
    {((3, 5, 0), (5, 8, 0)): 1},
    {((4, 7, 0), (7, 9, 1)): 1},
    {((5, 8, 0), (8, 7, 0)): 1},
    {((6, 7, 0), (7, 9, 0)): 1},
    {((8, 7, 0), (7, 9, 1)): 1, ((8, 7, 0), (7, 9, 0)): -1},
]


def parallelExample(swapKeys = False):
    """The algebra above, optionally with the parallel pair's names exchanged.

    Exchanging the two keys of the bundle `7 -> 9` is a renaming and nothing
    else: it is the same algebra, written down differently, which is exactly the
    ambiguity the procedure's build order creates between two routes to one
    quiver.
    """
    def rename(arrow):
        tail, head, key = arrow
        if swapKeys and (tail, head) == (7, 9):
            return (tail, head, 1 - key)
        return arrow

    quiver = qm.pathAlgebra.PathAlgebra().quiver.__class__()
    for tail, head, key in PARALLEL_ARROWS:
        quiver.add_edge(*rename((tail, head, key))[:2], key = rename((tail, head, key))[2])
    relations = [
        {tuple(rename(arrow) for arrow in path): coefficient
         for path, coefficient in relation.items()}
        for relation in PARALLEL_RELATIONS
    ]
    return procedure.toPathAlgebra(quiver, relations)


# -- the canonical key ----------------------------------------------------

def test_key_is_stable_and_separates():
    """Same algebra, same key; a changed relation, a changed key."""
    first = line_algebra(7, "22222")
    again = line_algebra(7, "22222")
    other = line_algebra(7, "22230")
    assert fingerprint.canonicalKey(first) == fingerprint.canonicalKey(again)
    assert fingerprint.canonicalKey(first) != fingerprint.canonicalKey(other)


def test_key_survives_a_deepcopy():
    """The search deepcopies at every step, so the key must not see identity."""
    alg = line_algebra(6, "3030")
    assert fingerprint.canonicalKey(alg) == fingerprint.canonicalKey(copy.deepcopy(alg))


def test_key_agrees_with_quiverKey_where_that_one_answers():
    """Where there are no parallel arrows the two must partition alike.

    `quiverKey` is the exact key the meeting-point machinery already trusts, so
    a canonical key that split or merged anything it does not would be a
    regression in the part of the search that is already right.
    """
    algebras = [line_algebra(7, name) for name in
                ["00000", "22222", "22230", "30300", "40000", "05000"]]
    plain = [search.quiverKey(alg) for alg in algebras]
    canonical = [fingerprint.canonicalKey(alg) for alg in algebras]
    assert all(key is not None for key in plain)
    for i in range(len(algebras)):
        for j in range(len(algebras)):
            assert (plain[i] == plain[j]) == (canonical[i] == canonical[j])


def test_parallel_arrow_renaming_does_not_change_the_key():
    """The case the old key had to refuse, and the reason this module exists."""
    straight, swapped = parallelExample(), parallelExample(swapKeys = True)
    assert arrowPaths.hasParallelArrows(straight.quiver)
    # The two really are written down differently: the raw arrow relations
    # differ, which is what made a search unable to meet at such a node.
    assert (procedure.relationsFrom(straight) != procedure.relationsFrom(swapped))
    assert search.quiverKey(straight) is None
    assert fingerprint.canonicalKey(straight) == fingerprint.canonicalKey(swapped)


def test_parallel_bundles_and_their_cost():
    """One bundle of two arrows, so two relabelings -- the case seen in practice."""
    straight = parallelExample()
    assert fingerprint.arrowBundles(straight.quiver) == {(7, 9): [0, 1]}
    assert fingerprint.relabelingCost(straight.quiver) == 2
    assert fingerprint.relabelingCost(line_algebra(7, "22222").quiver) == 1


def test_the_cap_refuses_rather_than_guesses():
    """Over the cap the key is None, and None must never read as a match."""
    straight = parallelExample()
    assert fingerprint.canonicalKey(straight, cap = 1) is None
    visited = fingerprint.Visited(cap = 1)
    assert not visited.seen(straight, 3)
    assert not visited.seen(straight, 3)     # refused twice, skipped neither
    assert visited.refused == 2
    assert visited.hits == 0


def test_a_relation_and_its_negation_are_one_relation():
    """`p - q = 0` and `q - p = 0` generate the same ideal, so one key."""
    positive = fingerprint._normaliseSign(((("a",), -1), (("b",), 1)))
    negative = fingerprint._normaliseSign(((("a",), 1), (("b",), -1)))
    assert positive == negative


# -- the digest -----------------------------------------------------------

def test_digest_is_deterministic_and_sized():
    key = fingerprint.canonicalKey(line_algebra(7, "22222"))
    assert fingerprint.digest(key) == fingerprint.digest(key)
    assert fingerprint.digest(key) < 2 ** 128
    assert fingerprint.digest(key, 64) < 2 ** 64
    assert fingerprint.digest(key) != fingerprint.digest(
        fingerprint.canonicalKey(line_algebra(7, "22230")))


def test_digest_passes_a_refusal_through():
    """A node the key refused must not acquire one by being hashed."""
    assert fingerprint.digest(None) is None


def test_digest_rejects_a_silly_width():
    with pytest.raises(ValueError):
        fingerprint.digest(fingerprint.canonicalKey(line_algebra(5, "000")), 17)


def test_collision_chance_is_the_birthday_bound():
    """Sanity, and the numbers the module docstring quotes."""
    assert fingerprint.collisionChance(1) == 0.0
    assert fingerprint.collisionChance(10 ** 8, 128) < 1e-22
    assert 1e-4 < fingerprint.collisionChance(10 ** 8, 64) < 1e-3
    assert (fingerprint.collisionChance(10 ** 6, 64)
            > fingerprint.collisionChance(10 ** 6, 128))


# -- the visited set ------------------------------------------------------

def test_visited_is_depth_aware():
    """A node met shallow must be walked again when met with more depth left.

    This is the whole correctness condition of the dedup: a plain visited set
    would skip the second visit and lose every branch below it.
    """
    alg = line_algebra(6, "3030")
    visited = fingerprint.Visited()
    assert not visited.seen(alg, 2)          # first sight, two levels left
    assert visited.seen(alg, 2)              # same depth, nothing new below
    assert visited.seen(alg, 1)              # shallower, already covered
    assert not visited.seen(alg, 5)          # deeper, three levels never walked
    assert visited.seen(alg, 5)


def test_visited_counts_add_up():
    alg = line_algebra(6, "3030")
    visited = fingerprint.Visited()
    visited.seen(alg, 3)
    visited.seen(alg, 3)
    summary = visited.summarise()
    assert summary['nodes'] == 2
    assert summary['distinct'] == 1
    assert summary['skipped'] == 1
    assert len(visited) == 1


def test_visited_can_store_digests_instead_of_keys():
    alg = line_algebra(6, "3030")
    visited = fingerprint.Visited(digestBits = 128)
    assert not visited.seen(alg, 3)
    assert visited.seen(alg, 3)
    assert visited.summarise()['digestBits'] == 128


# -- the deduped walk -----------------------------------------------------

DEDUPE_CASES = [
    (5, "000", 4), (5, "300", 4), (5, "230", 4),
    (6, "3030", 4), (6, "2300", 4), (6, "0400", 4),
    (7, "22222", 5), (7, "30300", 5), (7, "05000", 4),
]


@pytest.mark.parametrize("length,relLengths,depth", DEDUPE_CASES)
def test_dedup_reaches_exactly_what_the_plain_walk_reaches(length, relLengths, depth):
    """The routes are deduplicated; the destinations are not.

    Checked three ways at once, because the three are what callers read: the
    lines collected, the hereditary forms collected, and the set of distinct
    quivers the visitor is shown.
    """
    def walk(visited):
        alg = nk.LinearNakayamaAlgebra(length, [int(c) for c in relLengths])
        collected, hereditary, keys = [], [], set()
        search.mutationSearchDepthFirst(
            alg, depth, [], 'test', printOutput = False,
            collected = collected, collectedHereditary = hereditary,
            visitor = lambda pathAlg, mv: keys.add(fingerprint.canonicalKey(pathAlg)),
            visited = visited)
        lines = {(tuple(sorted(pathAlg.quiver.edges())),
                  tuple(sorted(tuple(sorted(tuple(p) for p in rel)) for rel in pathAlg.rels)))
                 for pathAlg, _path, _numbering in collected}
        return lines, {form for form, _quipu, _path in hereditary}, keys

    plain = walk(None)
    deduped = walk(fingerprint.Visited())
    assert deduped[0] == plain[0]
    assert deduped[1] == plain[1]
    assert deduped[2] == plain[2]


def test_dedup_actually_skips_something():
    """A dedup that reached the same answers by doing the same work would be
    correct and pointless, so pin that it is not that."""
    alg = line_algebra(7, "22222")
    visited = fingerprint.Visited()
    search.mutationSearchDepthFirst(alg, 5, [], 'test', printOutput = False,
                                    visited = visited)
    summary = visited.summarise()
    assert summary['skipped'] > 0
    assert summary['ratio'] > 1.0


def test_dedup_with_digests_reaches_the_same_answers():
    """At this size a 128-bit digest cannot collide, so it must agree exactly."""
    def walk(visited):
        alg = line_algebra(7, "22222")
        collected = []
        search.mutationSearchDepthFirst(alg, 5, [], 'test', printOutput = False,
                                        collected = collected, visited = visited)
        return {(tuple(sorted(p.quiver.edges())),
                 tuple(sorted(tuple(sorted(tuple(x) for x in r)) for r in p.rels)))
                for p, _a, _b in collected}

    assert walk(fingerprint.Visited(digestBits = 128)) == walk(None)


# -- the sign gauge -------------------------------------------------------

def gaugeExample(sign):
    """One algebra, written with the zero relation as `p = 0` or as `-p = 0`.

    The pair a depth-4 walk out of the n = 7 LNA `30300` reaches along two
    routes.  They are the same ideal -- a generator and its negative generate
    the same thing -- and the procedure nonetheless gives different answers from
    them, which is research F-050 and the reason the key quotients by the gauge.
    """
    quiver = qm.pathAlgebra.PathAlgebra().quiver.__class__()
    for tail, head in [(1, 2), (2, 6), (3, 4), (4, 1), (4, 5), (5, 6), (6, 7)]:
        quiver.add_edge(tail, head, key = 0)
    relations = [
        {((3, 4, 0), (4, 1, 0), (1, 2, 0)): sign},
        {((4, 1, 0), (1, 2, 0), (2, 6, 0)): 1, ((4, 5, 0), (5, 6, 0)): -1},
    ]
    return procedure.toPathAlgebra(quiver, relations)


def test_a_generator_and_its_negative_are_one_algebra():
    assert (fingerprint.canonicalKey(gaugeExample(1))
            == fingerprint.canonicalKey(gaugeExample(-1)))


def test_rescaling_an_arrow_does_not_change_the_key():
    """The two results the engine gives from those two presentations.

    They differ by the sign of one term, which is what rescaling an arrow does,
    so they are one algebra and must be one key.
    """
    from quivermutation import mutation, reduction
    results = [reduction.reducePathAlgebra(
                   mutation.quiverMutationAtVertex(gaugeExample(sign), 3))
               for sign in (1, -1)]
    # Different as written down ...
    assert (fingerprint.canonicalKey(results[0], gauge = False)
            != fingerprint.canonicalKey(results[1], gauge = False))
    # ... and one algebra once the gauge is taken out.
    assert fingerprint.canonicalKey(results[0]) == fingerprint.canonicalKey(results[1])


def test_the_gauge_does_not_merge_different_quivers():
    """Quotienting the signs must not start merging things that differ elsewhere."""
    def commuting(sign):
        quiver = qm.pathAlgebra.PathAlgebra().quiver.__class__()
        for tail, head in [(1, 2), (2, 4), (1, 3), (3, 4)]:
            quiver.add_edge(tail, head, key = 0)
        return procedure.toPathAlgebra(quiver, [
            {((1, 2, 0), (2, 4, 0)): 1, ((1, 3, 0), (3, 4, 0)): sign}])

    # The arrow 2 -> 4 lies in one path and not the other and nothing else
    # constrains it, so rescaling it carries one to the other: one algebra.
    assert (fingerprint.canonicalKey(commuting(1))
            == fingerprint.canonicalKey(commuting(-1)))
    # A different quiver is still a different key, gauge or no gauge.
    assert (fingerprint.canonicalKey(commuting(1))
            != fingerprint.canonicalKey(line_algebra(4, "20")))


def test_the_gauge_shrinks_the_walk_further():
    """Quotienting must leave strictly fewer distinct nodes on a walk that has
    the sign ambiguity in it at all."""
    alg = line_algebra(7, "30300")
    plain, quotiented = set(), set()
    search.mutationSearchDepthFirst(
        alg, 5, [], 'test', printOutput = False,
        visitor = lambda p, mv: (plain.add(fingerprint.canonicalKey(p, gauge = False)),
                                 quotiented.add(fingerprint.canonicalKey(p))))
    assert len(quotiented) < len(plain)


def test_the_coset_reduction_is_canonical_in_coset_and_minimal():
    """The three properties `_reduceModuloSpan` has to have.

    Canonical: two sign vectors reachable from one another reduce to the same
    thing, which is what makes the key well defined.  In the coset: the
    representative is a presentation the gauge can actually reach, not an
    invention.  Minimal: it is the least element of the coset, which is what
    makes it independent of the order the generators arrive in.

    Brute-forced against the whole coset, which is affordable at eight
    generators and is the only check that does not just restate the code.
    """
    import random
    rng = random.Random(99)
    for _ in range(400):
        generators = [rng.randrange(0, 1 << 12)
                      for _ in range(rng.randrange(1, 9))]
        vector = rng.randrange(0, 1 << 12)

        coset = {0}
        for generator in generators:
            coset |= {value ^ generator for value in coset}

        reduced = fingerprint._reduceModuloSpan(vector, generators)
        assert (reduced ^ vector) in coset
        assert reduced == min(vector ^ value for value in coset)

        elsewhere = vector
        for generator in generators:
            if rng.random() < 0.5:
                elsewhere ^= generator
        assert fingerprint._reduceModuloSpan(elsewhere, generators) == reduced
