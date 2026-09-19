"""Two families that could carry the classes the quipu theorem does not name.

The quipu theorem covers the LNAs whose relations are almost separate, and names
their classes by a tree with no relations.  Everything else has to be classified
some other way, and these are the two candidate families:

* every **tree**, not only the quipus -- `treeSearch`;
* every **quipu with relations**, in every orientation -- `quipuRelations`.

The tests pin what the searches found, and, more importantly, the checks that say
the searches are asking the right question: that the quipu theorem's own case
comes out of the enumeration unchanged, and that an algebra known to be in a
class by mutation is one the enumeration produces.
"""

import pytest

import quivermutation as qm
from quivermutation import coxeterTables as ct
from quivermutation import nakayama as nk
from quivermutation import quipuRelations as qr
from quivermutation import treeSearch as ts


# -- the status of every LNA ---------------------------------------------

@pytest.mark.parametrize("length, expected", [
    (6, {ct.QUIPU: 42, ct.NOT_QUIPU: 0, ct.UNPLACED: 0}),
    (7, {ct.QUIPU: 132, ct.NOT_QUIPU: 0, ct.UNPLACED: 0}),
    (8, {ct.QUIPU: 429, ct.NOT_QUIPU: 0, ct.UNPLACED: 0}),
    (9, {ct.QUIPU: 1421, ct.NOT_QUIPU: 9, ct.UNPLACED: 0}),
])
def test_how_many_lnas_the_known_moves_place_in_a_quipu_class(length, expected):
    """The counts of research F-032, recomputed through the polynomial tables."""
    assert ct.summary(length) == expected


def test_the_nine_lnas_outside_a_quipu_class_at_length_nine():
    """F-011's two classes: `3033030` and the eight-member `C(2,4,4)` orbit."""
    status = ct.lnaStatus(9)
    outside = sorted(ct.className(relLengths) for relLengths, value in status.items()
                     if value != ct.QUIPU)
    assert outside == ['3033030', '3345000', '3455000', '3505000', '4444400',
                       '4550400', '5040330', '5050030', '5504030']


# -- trees ----------------------------------------------------------------

@pytest.mark.parametrize("order, trees, notQuipus", [
    (9, 47, 29), (10, 106, 70), (11, 235, 171),
])
def test_how_many_trees_of_an_order_are_not_quipus(order, trees, notQuipus):
    rows = ts.treeReport(order)
    assert len(rows) == trees
    assert sum(1 for row in rows if not row['isQuipu']) == notQuipus


@pytest.mark.parametrize("order", [9, 10, 11])
def test_no_tree_outside_the_quipu_shape_matches_an_unclassified_lna(order):
    """The extension of F-031 to every degree: still nothing.

    A non-quipu tree that shared a Coxeter polynomial with an LNA in no quipu
    class would be a second hereditary family behind the classification.  There
    is none at these orders -- and the matches that do occur are all against LNAs
    the moves place in a quipu class, which the tree is cospectral with and so
    cannot be derived equivalent to.
    """
    assert ts.leads(order) == []


@pytest.mark.parametrize("order, cospectral", [(9, 3), (10, 0), (11, 7)])
def test_how_many_non_quipu_trees_are_cospectral_with_a_quipu(order, cospectral):
    assert len(ts.cospectralWithAQuipu(order)) == cospectral


# -- quipus with relations ------------------------------------------------

def test_the_line_is_a_quipu_and_is_left_out_of_the_enumeration():
    """One orientation of the cordless quipu is the line, whose algebras are LNAs."""
    order = 6
    linear = [orientation for _parameters, edges, orientation, _automorphisms
              in qr.orientedQuipus(order)
              if qr.isLinearlyOriented(order, edges, orientation)]
    assert linear == [(0, 0, 0, 0, 0)]


def test_a_quipu_algebra_reached_by_mutation_is_one_the_enumeration_produces():
    """The end-to-end check: mutation and enumeration must meet.

    One left-hand mutation of `3033030` -- the LNA of length 9 that is in no
    quipu class and is not piecewise hereditary -- gives a quipu quiver with
    relations, so it is in that LNA's class by construction.  Its certificate
    must therefore be one the enumeration produces, under exactly that LNA's
    Coxeter polynomial.
    """
    algebra = nk.LinearNakayamaAlgebra.fromClassName('3033030')
    mutated = qm.quiverMutationAtVertices(algebra, [2])
    assert qr.certificate(mutated) is not None
    result = qr.search(9, minArrows = 2, statuses = (ct.NOT_QUIPU,),
                       keepPerKey = None, dedupe = True)
    certificates = {match['certificate']: match for match in result['matches']}
    assert qr.certificate(mutated) in result['certificates']
    assert certificates[qr.certificate(mutated)]['lnas'] == (('3033030', ct.NOT_QUIPU),)


@pytest.mark.slow
def test_every_class_outside_the_quipu_theorem_at_length_nine_has_quipu_members():
    """Both classes the theorem misses at n = 9 contain quipus with relations.

    Walked from the other side: every quiver a mutation search out of the LNA
    reaches is in its class, so the quipu quivers among them are members with a
    path to prove it.  Every one of them must be in the enumeration, and the two
    classes must both have some.
    """
    result = qr.search(9, minArrows = 2, statuses = (ct.NOT_QUIPU, ct.UNPLACED),
                       keepPerKey = 0, dedupe = True)
    certificates = result['certificates']
    for className in ('3033030', '3345000'):
        algebra = nk.LinearNakayamaAlgebra.fromClassName(className)
        reached = qr.reachedQuipuAlgebras(algebra, 3)
        confirmed = {key for key in reached
                     if not _isLinearlyOrientedCertificate(key)}
        assert confirmed
        assert confirmed <= certificates


def _isLinearlyOrientedCertificate(certificate):
    """Whether a certificate's arrows are 0 -> 1 -> ... -> n - 1, i.e. an LNA."""
    arrows = certificate[0]
    return all(head == tail + 1 for tail, head in arrows)


@pytest.mark.parametrize("order, checked, agreeing", [(6, 543, 300), (7, 4160, 2138)])
def test_relations_of_two_arrows_are_not_free_on_a_quipu(order, checked, agreeing):
    """`corollary:lengthtworelations` does not extend past the line quiver.

    On an LNA a relation of two arrows can be deleted without leaving the derived
    equivalence class (F-028).  On a quipu with a branch it usually cannot: the
    Coxeter polynomial changes, and the polynomial is a derived invariant, so this
    is a refutation and not a gap in the evidence.
    """
    seen, kept, failures = qr.freeRelationCheck(order)
    assert (seen, kept) == (checked, agreeing)
    assert failures


def test_two_arrow_relations_are_free_on_the_line_itself():
    """The control: on the linearly oriented line the corollary does hold."""
    order = 7
    for parameters, edges, orientation, _automorphisms in qr.orientedQuipus(order):
        if not qr.isLinearlyOriented(order, edges, orientation):
            continue
        succ = qr.successors(order, edges, orientation)
        paths = qr.directedPaths(order, succ)
        baseMask, candidates, killMasks, comparable = qr.relationData(order, paths, 2)
        checked = 0
        for mask, chosen in qr.relationSets(baseMask, killMasks, comparable):
            shortOnes = [i for i in chosen if len(candidates[i]) - 1 == 2]
            if not shortOnes:
                continue
            stripped = baseMask
            for index in chosen:
                if index not in shortOnes:
                    stripped &= ~killMasks[index]
            from quivermutation import invariants as inv
            assert inv.coxeterCoefficients(qr.cartanFromMask(mask, order)) == \
                inv.coxeterCoefficients(qr.cartanFromMask(stripped, order))
            checked += 1
        assert checked > 50
        return
    raise AssertionError("the line was not in the enumeration")


# -- the command line -----------------------------------------------------

def test_the_three_reports_run(capsys):
    """A smoke test of `families.py`, at orders small enough to be instant."""
    import families

    assert families.main(["trees", "7"]) == 0
    assert "no non-quipu tree shares a polynomial" in capsys.readouterr().out

    assert families.main(["quipus", "7", "--min-arrows", "2"]) == 0
    printed = capsys.readouterr().out
    # Every LNA of length 7 is in a quipu class, so there is nothing to match.
    assert "matching against 0 LNAs" in printed

    assert families.main(["free", "4"]) == 0
    assert "4 do not" in capsys.readouterr().out


def test_the_members_report_prints_what_a_walk_proves(capsys):
    """`members` prints only algebras with a mutation path behind them."""
    import families

    assert families.main(["members", "9", "--depth", "2", "--show", "1"]) == 0
    printed = capsys.readouterr().out
    assert "9 LNAs in no quipu class" in printed
    assert "3033030:" in printed


# -- relation-free sightings ----------------------------------------------

def test_a_search_records_the_relation_free_quivers_it_reaches():
    """The sightings hook fires for every relation-free quiver, not just the first.

    `hereditaryFormsReachedFrom` keeps one entry per underlying graph; the
    sightings keep every visit, which is the point -- how many different ones a
    run sees, and whether any of them is not a tree, are questions the deduped
    answer cannot be asked.
    """
    from quivermutation import search as se

    algebra = nk.LinearNakayamaAlgebra.fromClassName('300')
    with se.relationFreeSightings() as sightings:
        forms = se.hereditaryFormsReachedFrom(algebra, 4)
    assert forms                                   # it does reach one
    assert len(sightings) > len(forms)             # and more than once
    summary = se.summariseSightings(sightings)
    assert summary['sightings'] == len(sightings)
    assert summary['distinct'] == len(forms)
    assert summary['quipus'] == len(sightings)     # D_5, every time
    assert summary['notTrees'] == 0
    assert summary['oddities'] == []
    for sighting in sightings:
        assert sighting['isTree'] and sighting['isQuipu']
        assert not sighting['hasOrientedCycle']
        assert sighting['parallelArrows'] == 0


def test_nothing_is_recorded_when_no_sink_is_open():
    """The hook costs nothing when it is not asked for."""
    from quivermutation import search as se

    assert se._SIGHTING_SINKS == []
    se.hereditaryFormsReachedFrom(nk.LinearNakayamaAlgebra.fromClassName('300'), 3)
    assert se._SIGHTING_SINKS == []


def test_the_sightings_are_written_as_json_lines(tmp_path, capsys):
    """`classify.py --sightings FILE` writes one JSON object per sighting."""
    import json

    import classify
    from quivermutation import search as se

    with se.relationFreeSightings() as sightings:
        se.hereditaryFormsReachedFrom(nk.LinearNakayamaAlgebra.fromClassName('300'), 3)
    target = tmp_path / "sightings.jsonl"
    counts = classify.write_sightings(sightings, str(target))
    written = [json.loads(line) for line in target.read_text().splitlines()]
    assert len(written) == len(sightings) == counts['sightings']
    assert written[0]['quipu'] == 'P^(2)_(1,1)'
    assert "relation-free quivers reached" in capsys.readouterr().out


# -- meeting in the middle, and orbits of one LNA --------------------------

def test_two_searches_that_meet_prove_a_mutation_equivalence():
    """A meeting point is a quiver both sides reach, so both are in one class.

    Vertex labels do not move under mutation, so the test is equality of
    labelled quivers rather than isomorphism, and a meeting is exact.
    """
    from quivermutation import search as se

    algebra = nk.LinearNakayamaAlgebra.fromClassName('300')
    mutated = qm.quiverMutationAtVertices(algebra, [4, 1])
    meetings = se.meetingPoints(algebra, mutated, 2)
    assert meetings
    key, fromFirst, fromSecond = meetings[0]
    assert key == se.quiverKey(qm.quiverMutationAtVertices(algebra, fromFirst))
    assert key == se.quiverKey(qm.quiverMutationAtVertices(mutated, fromSecond))


def test_meeting_in_the_middle_reaches_twice_as_far():
    """The pair the known moves miss at n = 8 meets at 3 mutations from each side.

    One-sided, that is a depth-6 search; the cost of a search is exponential in
    the depth, so this is the same reach for the square root of the work. F-041.
    """
    from quivermutation import search as se

    first = nk.LinearNakayamaAlgebra.fromClassName('230300')
    second = nk.LinearNakayamaAlgebra.fromClassName('030300')
    assert not se.meetingPoints(first, second, 2, alsoDual = True)
    meetings = se.meetingPoints(first, second, 3, alsoDual = True)
    assert meetings
    assert len(meetings[0][1]) + len(meetings[0][2]) == 6


def test_every_single_two_arrow_deletion_at_length_eight_is_a_mutation():
    """H-012's question, one relation at a time, settled at n = 8 (F-041).

    The whole strip is a composition of single deletions, so this says every LNA
    of length 8 is mutation equivalent to its stripped form.
    """
    from quivermutation import freeMoves as fm
    from quivermutation import lnaMoves as lm
    from quivermutation import search as se

    lnas, orbits = fm.derivedOrbits(8, rules = lm.ALL_MOVES, free = False,
                                    edges = True, doubles = True)
    where = {member: key for key, members in orbits.items() for member in members}
    open_ = []
    for lna in lnas:
        for start, arrows in enumerate(lna):
            if arrows != 2:
                continue
            without = list(lna)
            without[start] = 0
            without = tuple(without)
            if where[lna] != where[without]:
                open_.append((lna, without))
    assert len(open_) == 10                      # the known moves leave ten
    for lna, without in open_:
        assert se.meetingPoints(nk.LinearNakayamaAlgebra(8, list(lna)),
                                nk.LinearNakayamaAlgebra(8, list(without)),
                                3, alsoDual = True, limit = 1)


def test_the_orbit_of_one_lna_sits_inside_the_whole_length_partition():
    """`orbitOf` walks out of one row; `derivedOrbits` partitions every row.

    The walk is one-way, so it is contained in the class and need not exhaust it:
    `000000` has no relation for any move to act on, while plenty of rows
    double-mutate onto it.  That asymmetry is why the function is documented as a
    proof of reachability and not of its absence.
    """
    from quivermutation import freeMoves as fm

    lnas, orbits = fm.derivedOrbits(8, free = False, edges = True, doubles = True)
    byMember = {member: key for key, members in orbits.items() for member in members}
    for relLengths in [(3, 0, 3, 0, 3, 0), (2, 3, 0, 3, 0, 0), (0, 0, 0, 0, 0, 0)]:
        walked = fm.orbitOf(8, relLengths)
        assert walked <= set(orbits[byMember[relLengths]])
    # A row with relations does reach the whole of its class here.
    crowded = (2, 3, 0, 2, 0, 0)
    assert fm.orbitOf(8, crowded) <= set(orbits[byMember[crowded]])
    assert len(fm.orbitOf(8, crowded)) > 1


def test_the_barricade_shape_first_fits_at_length_thirteen():
    """Two heavy clusters walled around a two-arrow relation need 13 vertices.

    The configuration H-012's doubts are about, and F-040's point: it does not
    exist at any length this project has classified.
    """
    from quivermutation import freeMoves as fm
    from quivermutation import overlap as ov

    barricade = (3, 3, 0, 0, 0, 2, 0, 0, 3, 3, 0)
    algebra = nk.LinearNakayamaAlgebra(13, list(barricade))
    assert algebra.length == 13
    runs = [run for run in ov.overlapRuns(barricade, threshold = 2) if len(run) > 1]
    assert len(runs) == 2
    # and the known moves do get the two-arrow relation out of it
    assert fm.stripLengthTwo(barricade) in fm.orbitOf(13, barricade)


def test_meeting_in_the_middle_joins_a_barricade_a_one_way_walk_misses():
    """`movesJoin` walks from both ends, which is what the long quivers need.

    The barricades at `n = 15` and beyond have orbits of tens of thousands of
    rows, so a one-way walk under any affordable cap reports nothing and proves
    nothing.  Meeting in the middle finds the join immediately.  E-037.
    """
    from quivermutation import freeMoves as fm

    barricade = (3, 3, 0, 0, 0, 2, 0, 0, 0, 4, 4, 0, 0)
    strip = fm.stripLengthTwo(barricade)
    assert strip not in fm.orbitOf(15, barricade, limit = 400)
    meet = fm.movesJoin(15, barricade, strip, limit = 40000)
    assert meet is not None
    # the meeting row is genuinely reached from both ends
    assert meet in fm.orbitOf(15, barricade, limit = 40000, target = meet)
    assert meet in fm.orbitOf(15, strip, limit = 40000, target = meet)


def _reachesAlmostSeparate(length, relLengths, limit = 60000):
    """Does the full move table carry this LNA to one the quipu theorem names?"""
    from quivermutation import freeMoves as fm
    from quivermutation import lnaMoves as lm
    from quivermutation import overlap as ov

    seen, frontier = {tuple(relLengths)}, [tuple(relLengths)]
    while frontier and len(seen) < limit:
        current = frontier.pop()
        if ov.isAlmostSeparate(length, list(current)):
            return True
        for name in fm.movesFrom(length, current, lm.ALL_MOVES, False, True, True):
            if name not in seen:
                seen.add(name)
                frontier.append(name)
    return False


def test_where_a_core_sits_decides_it_and_not_how_much_it_overlaps():
    """The `45` core is reached from the ends and nowhere else.  F-042.

    A four-arrow relation and a five-arrow one sharing three arrows: one overlap,
    two relations, the same at every placement.  Against the source, or within one
    arrow of the sink, the moves carry it to an almost separate LNA; anywhere else
    they do not, however much free space is around it.
    """
    def place(length, offset):
        relLengths = [0] * (length - 2)
        relLengths[offset], relLengths[offset + 1] = 4, 5
        return tuple(relLengths)

    assert _reachesAlmostSeparate(9, place(9, 1))           # n = 9 has no room to escape
    for length in (11, 12):
        assert _reachesAlmostSeparate(length, place(length, 0))            # at the source
        assert _reachesAlmostSeparate(length, place(length, length - 7))   # at the sink
        assert not _reachesAlmostSeparate(length, place(length, 1))
        assert not _reachesAlmostSeparate(length, place(length, 2))
