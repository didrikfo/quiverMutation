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
