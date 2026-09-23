"""Placing an overlapping core in a long line, and the verdict it gets.

Two things are checked, and they are different in kind.  The first is that the
catalogue is made of real LNAs: `_rowFor` builds a row directly rather than by
filtering the enumerator, so it has to agree with the enumerator wherever the
enumerator can still be run, in both directions -- every row it builds is an
LNA, and every placement it refuses is genuinely not one.

The second is the answer itself.  F-042 published the `45` slide at lengths 9 to
13, so the task can be checked against a result that was arrived at
independently of it; a census instrument that cannot reproduce the one census
already done by hand is not worth running at a length where nothing is known.
"""

import argparse

import pytest

import batch
from quivermutation import nakayama as nk


WORDS = ['45', '54', '55', '33', '44', '504', '4004', '22', '26', '62', '345']


@pytest.mark.parametrize("length", [7, 8, 9, 10])
def test_every_row_the_catalogue_builds_is_an_LNA(length):
    every = set(nk.allRelationLengths(length))
    built = 0
    for word in batch._singleCores(3, 6) + WORDS:
        for offset in range(length):
            row = batch._rowFor(length, word, offset)
            if row is None:
                continue
            built += 1
            assert row in every, (word, offset, row)
    assert built > 0


@pytest.mark.parametrize("length", [8, 9])
def test_a_refused_placement_is_genuinely_not_an_LNA(length):
    """The other direction: `_rowFor` must not be refusing legal rows.

    Without this, a bug that refused everything awkward would pass the test
    above and quietly shrink the census to the easy placements.
    """
    every = set(nk.allRelationLengths(length))
    for word in WORDS:
        for offset in range(length):
            if batch._rowFor(length, word, offset) is not None:
                continue
            row = [0] * (length - 2)
            if offset < 0 or offset + len(word) > len(row):
                continue
            for position, letter in enumerate(word):
                if int(letter):
                    row[offset + position] = int(letter)
            assert tuple(row) not in every, (word, offset, tuple(row))


def test_the_catalogue_holds_only_heavily_overlapping_words():
    """A word the quipu theorem already names has nothing to ask, so it is out."""
    from quivermutation import overlap as ov

    for word in batch._singleCores(3, 6):
        code = [int(letter) for letter in word]
        assert ov.maxOverlap(code) >= 2, word
        assert code[0] and code[-1] and 1 not in code, word


#: F-042's table, read along each row: the `45` core at gap 0, 1, 2, ... from
#: the source.  Published 2026-09-18 from a hand slide, before this task existed.
FORTY_FIVE_SLIDE = {
    9:  ['inside', 'inside', 'inside'],
    10: ['inside', 'outside', 'inside', 'inside'],
    11: ['inside', 'outside', 'outside', 'inside', 'inside'],
    12: ['inside', 'outside', 'outside', 'outside', 'inside', 'inside'],
}


@pytest.mark.slow
@pytest.mark.parametrize("length", sorted(FORTY_FIVE_SLIDE))
def test_the_slide_reproduces_F042(length):
    verdicts = [batch._verdictFor(length, '45', offset,
                                  orbitLimit = 20000, joinLimit = 6000)['verdict']
                for offset in range(len(FORTY_FIVE_SLIDE[length]))]
    assert verdicts == FORTY_FIVE_SLIDE[length]


@pytest.mark.slow
@pytest.mark.parametrize("length", [11, 12])
def test_the_opposite_core_is_outside_at_the_mirrored_offsets(length):
    """F-042: reversing the arrows sends `45` at offset `o` to `504` at `n-7-o`.

    So the two slides are each other's mirror, and a task that got one right by
    accident would have to get the reflection right by the same accident.
    """
    fortyFive = {offset for offset in range(len(FORTY_FIVE_SLIDE[length]))
                 if FORTY_FIVE_SLIDE[length][offset] == 'outside'}
    fiveZeroFour = {offset for offset in range(length - 2)
                    if batch._rowFor(length, '504', offset) is not None
                    and batch._verdictFor(length, '504', offset, 20000, 6000)
                    ['verdict'] == 'outside'}
    assert fiveZeroFour == {length - 7 - offset for offset in fortyFive}


def test_a_verdict_is_never_outside_on_a_capped_orbit():
    """The E-037 guard: a walk that ran out of budget must not read as a result.

    Run with an orbit cap of 2, which every interesting row exceeds at once, and
    a join budget of 1, so neither step can conclude anything.  Nothing may come
    back `outside`: that verdict is reserved for an orbit that closed.
    """
    for offset in range(4):
        record = batch._verdictFor(12, '45', offset, orbitLimit = 2, joinLimit = 1)
        assert record['verdict'] != 'outside', record


def test_the_separated_targets_are_almost_separate_and_are_not_the_row():
    from quivermutation import overlap as ov

    row = batch._rowFor(15, '45', 3)
    targets = batch._separatedTargets(15, row)
    assert targets
    for target in targets:
        assert target != row
        assert ov.isAlmostSeparate(15, target)
        assert target in set(nk.allRelationLengths(15)) or len(target) == 13


def test_the_ledger_name_carries_every_parameter_that_decides_a_verdict():
    """Two runs that would answer differently must not share a ledger."""
    task = batch.CoresTask()
    parser = argparse.ArgumentParser()
    task.addArguments(parser)
    base = task.ledgerPath(parser.parse_args(["15"]))
    for flags in (["15", "--orbit-limit", "500"], ["15", "--join-limit", "50"],
                  ["15", "--max-word", "4"], ["15", "--gaps", ""],
                  ["14"], ["15", "--max-arrows", "5"],
                  ["15", "--pair-word", "3"]):
        assert task.ledgerPath(parser.parse_args(flags)) != base, flags


def test_no_two_parameter_sets_that_mean_different_work_share_a_ledger():
    """The general form of the test above, over the whole option surface.

    `--pair-word` was left out of the ledger name when this task was written,
    so `--max-word 4` and `--max-word 4 --pair-word 3` -- 2591 placements and
    17293 -- would have appended to the same file and each read the other's
    rows as its own finished work.
    """
    task = batch.CoresTask()
    parser = argparse.ArgumentParser()
    task.addArguments(parser)
    settings = [["13"], ["13", "--max-word", "4"], ["13", "--pair-word", "3"],
                ["13", "--max-arrows", "5"], ["13", "--gaps", "1"],
                ["13", "--gaps", ""], ["13", "--orbit-limit", "99"],
                ["13", "--join-limit", "99"], ["14"]]
    paths = {}
    for flags in settings:
        args = parser.parse_args(flags)
        path = task.ledgerPath(args)
        units = tuple(task.units(args))
        if path in paths:
            assert paths[path] == units, (flags, path)
        paths[path] = units
    assert len(paths) == len(settings)


def test_the_unit_list_is_stable_and_names_a_placement():
    task = batch.CoresTask()
    parser = argparse.ArgumentParser()
    task.addArguments(parser)
    args = parser.parse_args(["12", "--gaps", "1"])
    units = task.units(args)
    assert units == task.units(args)
    assert len(units) == len(set(units))
    for unit in units[:20]:
        word, offset = unit.split("@")
        assert batch._rowFor(12, word, int(offset)) is not None


# -- narrowing a census without splitting its ledger -----------------------
#
# The second overnight run left n = 13 at 62 cores and n = 15 at 125, because
# each length walked the catalogue at its own speed and the budget cut it
# wherever it had got to.  Only the overlap could be read, and reading a census
# against another length is the whole point of running one.  `--core-limit` and
# `--cores` make the cut deliberate instead; the ledger has to stay shared, or
# the second night would redo the first night's work.

def _coresArgs(*flags):
    task = batch.CoresTask()
    parser = argparse.ArgumentParser()
    task.addArguments(parser)
    return task, parser.parse_args(list(flags))


def test_core_limit_takes_a_prefix_of_the_catalogue_and_nothing_else():
    task, wide = _coresArgs("13", "--max-word", "3")
    _task, narrow = _coresArgs("13", "--max-word", "3", "--core-limit", "8")
    allUnits, someUnits = task.units(wide), task.units(narrow)
    assert someUnits == allUnits[:len(someUnits)], "a prefix and not a sample"
    assert 0 < len(someUnits) < len(allUnits)
    # Some of the first eight words may not fit the length at any offset, so
    # what is pinned is the prefix, not a count of eight.
    assert len({unit.split("@")[0] for unit in someUnits}) <= 8


def test_the_same_core_limit_asks_about_the_same_cores_at_every_length():
    # This is what makes two part-finished censuses comparable at all.
    words = []
    for length in (11, 13, 15, 17):
        task, args = _coresArgs(str(length), "--max-word", "3", "--core-limit", "12")
        words.append({unit.split("@")[0] for unit in task.units(args)})
    assert all(other == words[0] for other in words[1:])


def test_naming_cores_selects_exactly_those_and_keeps_catalogue_order():
    # Without the mirror: `504` is `45` reflected, so a mirrored run asks only
    # the `45` half of the pair, which is `test_reduced_walk`'s business.
    task, args = _coresArgs("13", "--max-word", "3", "--cores", "45,504,33",
                            "--no-mirror")
    units = task.units(args)
    assert {unit.split("@")[0] for unit in units} == {"45", "504", "33"}
    _task, wide = _coresArgs("13", "--max-word", "3", "--no-mirror")
    order = [unit for unit in task.units(wide) if unit in set(units)]
    assert units == order


def test_a_core_that_is_not_in_the_catalogue_is_refused_rather_than_ignored():
    task, args = _coresArgs("13", "--max-word", "3", "--cores", "45,9999")
    with pytest.raises(ValueError) as raised:
        task.units(args)
    assert "9999" in str(raised.value)


def test_narrowing_a_run_does_not_move_its_ledger():
    # A placement means the same thing however the run was narrowed, so a night
    # that does part of a census and a night that does the rest must write the
    # same file -- otherwise the second redoes the first.
    task, wide = _coresArgs("15", "--max-word", "4")
    _t, limited = _coresArgs("15", "--max-word", "4", "--core-limit", "40")
    _u, named = _coresArgs("15", "--max-word", "4", "--cores", "45")
    assert task.ledgerPath(wide) == task.ledgerPath(limited) == task.ledgerPath(named)


# -- reading a slide -------------------------------------------------------

def test_the_head_and_tail_are_the_inside_offsets_at_each_end():
    assert batch._headAndTail("iooooooii") == (1, 2, "all outside")
    assert batch._headAndTail("ooooooo") == (0, 0, "all outside")
    assert batch._headAndTail("iiiii") == (5, 0, "all outside")


def test_an_unfinished_or_undecided_interior_is_not_reported_as_settled():
    # The one reading that must never be produced by a half-done run.
    assert batch._headAndTail("ii..oooo")[2] == "not finished"
    assert batch._headAndTail("ooooo?ii")[2] == "not finished"
    assert batch._headAndTail("iooioooi")[2] == "HOLDS AN INSIDE"


# -- the verdict engine after the walk learned to stop early ---------------

def test_a_verdict_records_why_its_walk_stopped():
    record = batch._verdictFor(11, "45", 0, orbitLimit = 20000, joinLimit = 2000)
    assert record['verdict'] == 'inside'
    assert record['stopped'] == 'found'
    assert record['orbit'] <= 20000


def test_an_outside_verdict_still_requires_a_closed_orbit():
    record = batch._verdictFor(11, "45", 2, orbitLimit = 20000, joinLimit = 2000)
    assert record['verdict'] == 'outside'
    assert record['stopped'] == 'closed', \
        "outside may only be said of an orbit that emptied its frontier"


@pytest.mark.slow
def test_stopping_early_gives_the_verdicts_the_exhaustive_scan_gave():
    # Pins the economy to the answers F-042 published, at the length where the
    # `45` slide first has an outside in it.
    from quivermutation import freeMoves as fm, overlap as ov
    for offset in range(8):
        row = batch._rowFor(10, "45", offset)
        if row is None:
            continue
        whole = fm.orbitOf(10, row, free = True, edges = True, doubles = True,
                           limit = 20000)
        byScan = any(ov.isAlmostSeparate(10, member) for member in whole)
        record = batch._verdictFor(10, "45", offset, orbitLimit = 20000,
                                   joinLimit = 2000)
        assert (record['verdict'] == 'inside') == byScan, (offset, record)


def test_the_sampler_s_orbit_limit_is_in_its_ledger_name():
    # Same class of mistake as the one `--pair-word` made here: the limit
    # decides whether a draw is recorded as a leftover, so a ledger written
    # under one limit must not be read, or resumed, as though it answered
    # under another.
    task = batch.SampleTask()
    parser = argparse.ArgumentParser()
    task.addArguments(parser)
    tight = parser.parse_args(["15", "--orbit-limit", "20000"])
    loose = parser.parse_args(["15", "--orbit-limit", "200000"])
    assert task.ledgerPath(tight) != task.ledgerPath(loose)


def test_min_word_is_a_filter_that_leaves_out_the_narrower_catalogue():
    parser = argparse.ArgumentParser()
    task = batch.CoresTask()
    task.addArguments(parser)
    wide = parser.parse_args(["14", "--max-word", "5", "--gaps", ""])
    new = parser.parse_args(["14", "--max-word", "5", "--gaps", "", "--min-word", "5"])
    narrow = parser.parse_args(["14", "--max-word", "4", "--gaps", ""])
    assert task.ledgerPath(new) == task.ledgerPath(wide)
    # Together they ask everything the wide census asks.  They may ask a little
    # more: a five-letter word whose mirror is a four-letter word is asked here
    # because its mirror is not in this filtered list, and another ledger of
    # the length then answers it for nothing (E-051).
    assert set(task.units(new)) | set(task.units(narrow)) >= set(task.units(wide))
    assert all(len(unit.split("@")[0]) >= 5 for unit in task.units(new))
