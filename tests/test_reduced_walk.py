"""The free move walked both ways: an LNA and its stripped form are one state.

`corollary:lengthtworelations` makes a relation of two arrows free -- it can be
deleted or added -- and the plain walk (`free = True`) only ever deleted.  That
made the walk give two rows of one derived class different verdicts (E-049),
and made the core census ask every such pair twice.
"""

import argparse

import pytest

import batch
from quivermutation import freeMoves as fm
from quivermutation import nakayama as nk


@pytest.mark.parametrize("length", [7, 8])
def test_a_reduced_step_only_ever_names_reduced_rows(length):
    for row in list(nk.allRelationLengths(length))[:300]:
        for reached in fm.reducedMovesFrom(length, row):
            assert fm.isReduced(reached), (row, reached)
            assert reached in set(nk.allRelationLengths(length))


@pytest.mark.parametrize("length", [7, 8, 9])
def test_the_addable_vertices_are_exactly_the_admissible_additions(length):
    every = set(nk.allRelationLengths(length))
    for row in list(every)[:400]:
        for position, arrows in enumerate(row):
            if arrows:
                continue
            candidate = list(row)
            candidate[position] = 2
            assert ((position + 1) in fm.addableLengthTwo(length, row)) == \
                (tuple(candidate) in every), (row, position)


def test_an_lna_and_its_stripped_form_are_one_state():
    row = (2, 4, 5, 0, 0, 0, 0, 0, 0)
    stripped = fm.stripLengthTwo(row)
    first = fm.orbitReport(11, row, free = fm.REDUCED, limit = 5000)
    second = fm.orbitReport(11, stripped, free = fm.REDUCED, limit = 5000)
    assert first.rows == second.rows


def test_adding_a_two_arrow_relation_places_what_deleting_never_could():
    # E-049's smallest case: `404` at offset 1 of n = 11 has no move at all
    # under the plain walk, while `2404` at offset 0 -- the same LNA with a
    # two-arrow relation at the source -- reaches an almost separate row.
    plain = batch._verdictFor(11, "404", 1, 20000, 2000, free = True)
    alias = batch._verdictFor(11, "2404", 0, 20000, 2000, free = True)
    reduced = batch._verdictFor(11, "404", 1, 20000, 2000, free = fm.REDUCED)
    assert plain['verdict'] == 'outside' and plain['orbit'] == 1
    assert alias['verdict'] == 'inside'
    assert reduced['verdict'] == 'inside'


def _args(*flags):
    task = batch.CoresTask()
    parser = argparse.ArgumentParser()
    task.addArguments(parser)
    return task, parser.parse_args(list(flags))


def test_the_reduced_catalogue_holds_no_two_arrow_relation():
    task, args = _args("14", "--max-word", "4", "--walk", "reduced")
    for unit in task.units(args):
        word, offset = unit.split("@")
        assert fm.isReduced(batch._rowFor(14, word, int(offset))), unit


@pytest.mark.parametrize("length", [11, 14])
def test_every_alias_the_plain_catalogue_asks_is_a_reduced_placement(length):
    # The census loses nothing by dropping them: each is its stripped row, and
    # that row is already a unit of the reduced catalogue.
    task, plain = _args(str(length), "--max-word", "4", "--walk", "plain",
                        "--no-mirror")
    _t, reduced = _args(str(length), "--max-word", "4", "--walk", "reduced",
                        "--no-mirror")
    rows = {batch._rowFor(length, *_split(unit)) for unit in task.units(reduced)}
    aliases = 0
    for unit in task.units(plain):
        row = batch._rowFor(length, *_split(unit))
        if not fm.isReduced(row):
            aliases += 1
            assert fm.stripLengthTwo(row) in rows, unit
    assert aliases > 0


def _split(unit):
    word, offset = unit.split("@")
    return word, int(offset)


def test_asking_for_an_alias_under_the_reduced_walk_says_what_it_is():
    task, args = _args("13", "--max-word", "4", "--cores", "45,245",
                       "--walk", "reduced")
    with pytest.raises(ValueError) as raised:
        task.units(args)
    assert "245" in str(raised.value) and "two arrows" in str(raised.value)
    task, plain = _args("13", "--max-word", "4", "--cores", "45,245",
                        "--walk", "plain")
    assert any(unit.startswith("245@") for unit in task.units(plain))


def test_the_walk_is_in_both_ledger_names_and_plain_keeps_the_old_name():
    task, reduced = _args("15", "--max-word", "4", "--walk", "reduced")
    _t, plain = _args("15", "--max-word", "4", "--walk", "plain")
    assert task.ledgerPath(reduced) != task.ledgerPath(plain)
    assert task.ledgerPath(plain) == "logs/cores-n15-w4p2a6g123-o20000j6000.jsonl"
    sample = batch.SampleTask()
    parser = argparse.ArgumentParser()
    sample.addArguments(parser)
    assert sample.ledgerPath(parser.parse_args(["15", "--walk", "reduced"])) != \
        sample.ledgerPath(parser.parse_args(["15"]))
    assert sample.ledgerPath(parser.parse_args(["15"])) == \
        "logs/sample-n15-s0-d0-o20000.jsonl"


# -- the relation dual ------------------------------------------------------

def test_the_mirror_of_45_is_504_at_the_reflected_offset():
    for length in (11, 14, 17):
        for offset in range(length - 6):
            row = batch._rowFor(length, "45", offset)
            if row is None:
                continue
            assert batch._placementOf(batch._mirror(length, row)) == \
                ("504", length - 7 - offset)
            assert batch._mirror(length, batch._mirror(length, row)) == row


@pytest.mark.parametrize("walk", ["plain", "reduced"])
def test_a_mirrored_census_asks_one_row_of_each_dual_pair_and_loses_none(walk):
    task, whole = _args("13", "--max-word", "4", "--walk", walk, "--no-mirror")
    _t, halved = _args("13", "--max-word", "4", "--walk", walk)
    every = {batch._rowFor(13, *_split(unit)) for unit in task.units(whole)}
    asked = {batch._rowFor(13, *_split(unit)) for unit in task.units(halved)}
    assert asked < every
    for row in every - asked:
        assert batch._mirror(13, row) in asked, row
    # Filtering is not a different instrument: the ledger is the same one.
    assert task.ledgerPath(whole) == task.ledgerPath(halved)


def test_the_summary_reads_the_unasked_half_of_a_slide_in_the_mirror():
    import io
    task, args = _args("12", "--max-word", "4", "--cores", "45,504")
    verdicts = {0: 'inside', 1: 'outside', 2: 'outside', 3: 'outside',
                4: 'inside', 5: 'inside'}
    records = [{'result': {'core': '45', 'offset': offset, 'verdict': verdict,
                           'name': "".join(map(str, batch._rowFor(12, '45', offset)))}}
               for offset, verdict in verdicts.items()]
    out = io.StringIO()
    task.summarise(records, args, out)
    assert "504" in out.getvalue() and "iioooi" in out.getvalue()


# -- the shared walk ---------------------------------------------------------

def test_the_shared_walk_gives_the_reduced_verdict_and_spreads_it_to_the_class():
    walker = fm.SharedWalk(11)
    bare = batch._rowFor(11, "404", 1)
    alias = batch._rowFor(11, "2404", 0)
    verdict, how, _walked, _found, _promotes = walker.verdict(bare)
    assert verdict == 'inside' and how in ('plain', 'reduced')
    # The alias and the mirror are the same class, so neither is walked again.
    assert walker.verdict(alias)[:2] == ('inside', 'shared')
    assert walker.verdict(fm.mirrorRow(11, bare))[:2] == ('inside', 'shared')


def test_an_outside_verdict_is_a_closed_orbit_and_is_cached_row_by_row():
    walker = fm.SharedWalk(11)
    first = walker.verdict(batch._rowFor(11, "45", 1))
    assert first[:2] == ('outside', 'closed')
    # `45` at 2 is in the same closed orbit (F-051): its walk expands nothing new.
    second = walker.verdict(batch._rowFor(11, "45", 2))
    assert second[:2] == ('outside', 'closed')
    assert second[2] < first[2]


def test_a_later_certificate_promotes_an_earlier_outside_in_its_class():
    walker = fm.SharedWalk(11)
    row = batch._rowFor(11, "45", 1)
    assert walker.verdict(row, label = "45@1")[0] == 'outside'
    # Pretend a join found a certificate for a row of that class.
    assert walker.markInside(fm.mirrorRow(11, row)) == ["45@1"]
    assert walker.isInside(row)


def test_the_summary_applies_promotions_recorded_by_later_units():
    records = [
        {'unit': 'a', 'result': {'label': 'a', 'verdict': 'outside'}},
        {'unit': 'b', 'result': {'label': 'b', 'verdict': 'inside', 'promotes': ['a']}},
    ]
    latest = {result['label']: result['verdict']
              for result in batch._latest(records, ('verdict', 'inside'))}
    assert latest == {'a': 'inside', 'b': 'inside'}


@pytest.mark.slow
def test_a_shared_census_agrees_with_the_reduced_walk_placement_by_placement():
    # E-050 at n = 11, in miniature: every placement of the `--max-word 3`
    # catalogue, one process, against the reduced walk run alone on each.
    task, args = _args("11", "--max-word", "3", "--gaps", "")
    rows = [batch._rowFor(11, *_split(unit)) for unit in task.units(args)]
    walker = fm.SharedWalk(11)
    shared = [walker.verdict(row)[0] for row in rows]
    shared = ['inside' if walker.isInside(row) else verdict
              for row, verdict in zip(rows, shared)]
    alone = [batch._verdictFor(11, *_split(unit), 20000, 2000,
                               free = fm.REDUCED)['verdict']
             for unit in task.units(args)]
    # The reduced walk alone is asked of the unit's own row; an alias there can
    # lose to its stripped row only under the plain walk, never the reduced one.
    assert shared == alone
