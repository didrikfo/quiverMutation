"""Drawing LNAs uniformly at lengths too long to enumerate.

The claim is that the draw is *uniform over the same set the enumerator
produces*, so both halves are checked against `lines`: the count against how
many it produces, and the support against which ones.  A sampler that drew only
the sparse LNAs would pass a count test and fail this.
"""

import collections
import random

import pytest

from quivermutation import lines
from quivermutation import sampling


def everyLNA(length):
    """The enumerator's answer, as relation-length tuples."""
    return {tuple(lines.relationStringToLineRelLengths(length, lines.relSetToString(relSet)))
            for relSet in lines.generateAllPossibleLineRelations(length)}


CATALAN = [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796, 58786]


@pytest.mark.parametrize("length", range(1, 12))
def test_the_count_is_the_enumerator_s(length):
    assert sampling.countLNAs(length) == len(lines.generateAllPossibleLineRelations(length))


@pytest.mark.parametrize("length", range(1, 12))
def test_the_count_is_catalan(length):
    """Catalan(length - 1), which is what NOTES.md and the papers say."""
    assert sampling.countLNAs(length) == CATALAN[length - 1]


def test_the_count_reaches_lengths_the_enumerator_cannot():
    """The whole point: a number for a length nobody can list."""
    assert sampling.countLNAs(20) == 1767263190
    assert sampling.countLNAs(30) > 10 ** 15


@pytest.mark.parametrize("length", [5, 6, 7])
def test_every_LNA_is_drawn_and_none_that_is_not(length):
    """The support of the sampler is exactly the set of LNAs of the length."""
    rng = random.Random(20260919)
    counts = sampling.relationSetCounts(length)
    drawn = {tuple(sampling.sampleRelLengths(length, rng, counts))
             for _ in range(4000)}
    assert drawn == everyLNA(length)


@pytest.mark.parametrize("length", [5, 6])
def test_the_draw_is_uniform(length):
    """Every LNA within a fifth of its expected share over 40000 draws.

    A loose band on purpose: this is a test that the sampler is not *biased*,
    not a test of the random number generator, and the seed is fixed so it
    cannot start failing on a Tuesday.
    """
    rng = random.Random(4)
    counts = sampling.relationSetCounts(length)
    draws = 40000
    tally = collections.Counter(tuple(sampling.sampleRelLengths(length, rng, counts))
                                for _ in range(draws))
    expected = draws / sampling.countLNAs(length)
    assert min(tally.values()) > 0.8 * expected
    assert max(tally.values()) < 1.2 * expected


def test_a_drawn_LNA_is_a_legal_one():
    """Relations must fit: a relation at vertex v of k arrows needs v + k <= n."""
    from quivermutation import nakayama as nk
    rng = random.Random(11)
    for length in (12, 16, 20):
        counts = sampling.relationSetCounts(length)
        for _ in range(50):
            relLengths = sampling.sampleRelLengths(length, rng, counts)
            assert len(relLengths) == length - 2
            # Constructing it is the real test: the class raises on a relation
            # that does not fit.
            nk.LinearNakayamaAlgebra(length, list(relLengths))


def test_relations_have_increasing_starts_and_ends():
    """The structure the enumerator's recursion produces, and the DP counts."""
    rng = random.Random(5)
    for length in (9, 14, 18):
        counts = sampling.relationSetCounts(length)
        for _ in range(60):
            relations = sampling.sampleRelations(length, rng, counts)
            starts = [start for start, _end in relations]
            ends = [end for _start, end in relations]
            assert starts == sorted(set(starts))
            assert ends == sorted(set(ends))
            assert all(end - start >= 2 for start, end in relations)
            assert all(end <= length for _start, end in relations)


def test_a_draw_is_reproducible_and_independent_of_its_neighbours():
    """A resumed run must redraw exactly what the first run drew.

    Indexing the draws rather than streaming them is what lets a run be extended
    from 1000 to 2000 without redoing or disturbing the first 1000.
    """
    first = [sampling.drawFor(15, 7, index) for index in range(20)]
    again = [sampling.drawFor(15, 7, index) for index in range(20)]
    assert first == again
    # Drawing them in a different order changes nothing.
    assert sampling.drawFor(15, 7, 19) == first[19]
    # A different seed is a different sample.
    assert [sampling.drawFor(15, 8, index) for index in range(20)] != first


def test_draws_are_not_all_the_same():
    """A cheap guard against an index that is silently ignored."""
    drawn = {tuple(sampling.drawFor(16, 3, index)) for index in range(40)}
    assert len(drawn) > 30


# -- the probe ------------------------------------------------------------

def test_probe_settles_an_almost_separate_LNA_by_the_theorem():
    record = sampling.probe(7, [2, 0, 2, 0, 0])
    assert record['settledBy'] == 'theorem'
    assert record['name'] == '20200'
    assert record['orbit'] == 1


def test_probe_reports_the_overlap_profile():
    record = sampling.probe(9, [3, 5, 0, 5, 0, 0, 0])
    assert record['length'] == 9
    assert record['relations'] == 3
    assert record['maxOverlap'] >= 2
    assert record['settledBy'] in ('theorem', 'moves', 'leftover')


def test_probe_finds_the_known_n9_leftovers():
    """The two orbits `merges.py` is aimed at must come back as leftovers.

    This is the probe's calibration: it is the cheap pipeline asked one row at a
    time instead of a length at a time, and on a length where the exhaustive
    answer is known it has to agree with it.
    """
    for relLengths in [(3, 5, 0, 5, 0, 0, 0), (3, 0, 3, 3, 0, 3, 0)]:
        assert sampling.probe(9, relLengths)['settledBy'] == 'leftover'


def test_probe_is_json_serialisable():
    """It goes straight into a ledger line, so it must survive the round trip."""
    import json
    record = sampling.probe(9, [3, 0, 3, 3, 0, 3, 0])
    assert json.loads(json.dumps(record)) == record


# -- which kind of leftover it is ------------------------------------------
#
# `settledBy == 'leftover'` covers two different facts: an orbit that emptied
# its frontier without holding an almost separate row, and a walk that ran out
# of budget.  At n = 15 one leftover in six was the second kind, so a leftover
# rate read without the split is part rate and part cap.  E-037 is what reading
# one as the other costs.

def test_a_probe_says_whether_its_orbit_closed():
    record = sampling.probe(9, (0, 4, 5, 0, 0, 0, 0))
    assert record['orbitClosed'] in (True, False)


def test_a_leftover_found_under_a_tiny_limit_is_not_reported_as_closed():
    record = sampling.probe(13, (0, 0, 0, 4, 5, 0, 0, 0, 0, 0, 0), orbitLimit = 30)
    if record['settledBy'] == 'leftover':
        assert record['orbitClosed'] is False


def test_a_row_the_theorem_names_needs_no_walk_and_says_its_orbit_closed():
    record = sampling.probe(9, (2, 0, 0, 2, 0, 0, 0))
    assert record['settledBy'] == 'theorem'
    assert record['orbitClosed'] is True
    assert record['orbit'] == 1
