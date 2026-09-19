"""The ledger and the driver a long run is built on.

What has to hold for a run to be leaveable: a finished unit is on disk before
the next one starts, a resumed run does not redo it, a budget stops between
units rather than in the middle of one, and a ledger damaged by a kill is read
rather than refused.
"""

import json
import os
import sys

import pytest

from quivermutation import jobs


class Counting(jobs.Task):
    """A task whose unit is a number and whose work is squaring it."""

    name = 'counting'
    help = "squares numbers"

    def __init__(self, path, units = 5, boom = None):
        self.path = path
        self.count = units
        self.boom = boom

    def ledgerPath(self, args):
        return self.path

    def units(self, args):
        return [str(index) for index in range(self.count)]

    def run(self, unit, args):
        if self.boom is not None and unit == self.boom:
            raise RuntimeError("unit {0} was asked to fail".format(unit))
        return {'square': int(unit) ** 2}

    def summarise(self, records, args, out = None):
        out = sys.stdout if out is None else out
        print("squares: {0}".format(
            sorted(record['result']['square'] for record in records)), file = out)


class Slow(Counting):
    """The same, but each unit takes long enough for a budget to expire."""

    name = 'slow'

    def run(self, unit, args):
        import time
        time.sleep(0.25)
        return {'square': int(unit) ** 2}


@pytest.fixture
def ledgerPath(tmp_path):
    return str(tmp_path / "nested" / "ledger.jsonl")


# -- the ledger -----------------------------------------------------------

def test_a_missing_ledger_reads_as_empty(ledgerPath):
    ledger = jobs.Ledger(ledgerPath)
    assert ledger.records() == []
    assert ledger.done() == set()


def test_the_ledger_makes_its_directory(ledgerPath):
    jobs.Ledger(ledgerPath)
    assert os.path.isdir(os.path.dirname(ledgerPath))


def test_a_record_is_on_disk_before_the_next_one(ledgerPath):
    """The whole point of appending and flushing: a kill loses only what was in
    flight, so a reader opening the file mid-run sees the finished units."""
    ledger = jobs.Ledger(ledgerPath)
    ledger.record("a", {'value': 1})
    assert jobs.Ledger(ledgerPath).done() == {"a"}
    ledger.record("b", {'value': 2})
    assert jobs.Ledger(ledgerPath).done() == {"a", "b"}
    ledger.close()


def test_a_record_carries_its_result_and_a_timestamp(ledgerPath):
    with jobs.Ledger(ledgerPath) as ledger:
        ledger.record("a", {'value': 1})
    record = jobs.Ledger(ledgerPath).records()[0]
    assert record['unit'] == "a"
    assert record['result'] == {'value': 1}
    assert record['at']


def test_a_half_written_line_is_skipped_not_raised(ledgerPath):
    """A process killed mid-append leaves one broken line.  Refusing to start
    because of it would make a crash cost the whole run."""
    with jobs.Ledger(ledgerPath) as ledger:
        ledger.record("a", {'value': 1})
        ledger.record("b", {'value': 2})
    with open(ledgerPath, "a") as handle:
        handle.write('{"unit": "c", "resu')
    assert jobs.Ledger(ledgerPath).done() == {"a", "b"}


def test_blank_lines_are_skipped(ledgerPath):
    jobs.Ledger(ledgerPath)          # makes the directory
    with open(ledgerPath, "w") as handle:
        handle.write("\n" + json.dumps({'unit': 'a', 'result': {}}) + "\n\n")
    assert jobs.Ledger(ledgerPath).done() == {"a"}


# -- the driver -----------------------------------------------------------

def test_a_run_does_every_unit_and_records_it(ledgerPath, capsys):
    task = Counting(ledgerPath, units = 5)
    assert jobs.runTask(task, None) == 0
    assert jobs.Ledger(ledgerPath).done() == {"0", "1", "2", "3", "4"}


def test_a_unit_records_how_long_it_took(ledgerPath):
    jobs.runTask(Counting(ledgerPath, units = 2), None)
    for record in jobs.Ledger(ledgerPath).records():
        assert 'seconds' in record['result']


def test_a_rerun_does_nothing(ledgerPath, capsys):
    task = Counting(ledgerPath, units = 4)
    jobs.runTask(task, None)
    capsys.readouterr()
    assert jobs.runTask(task, None) == 0
    assert "nothing to do" in capsys.readouterr().out
    assert len(jobs.Ledger(ledgerPath).records()) == 4


def test_a_widened_run_does_only_the_new_units(ledgerPath, capsys):
    """Asking for more units must not redo the ones already in the ledger."""
    jobs.runTask(Counting(ledgerPath, units = 3), None)
    capsys.readouterr()
    jobs.runTask(Counting(ledgerPath, units = 6), None)
    output = capsys.readouterr().out
    assert "3 units already done" in output
    units = [record['unit'] for record in jobs.Ledger(ledgerPath).records()]
    assert sorted(units) == ["0", "1", "2", "3", "4", "5"]
    assert len(units) == len(set(units))


def test_a_budget_stops_between_units_and_exits_2(ledgerPath):
    """Exit 2 is the convention `overnight.py` restarts on, and what is finished
    must be on disk when it happens."""
    task = Slow(ledgerPath, units = 20)
    code = jobs.runTask(task, None, budgetHours = 0.5 / 3600)
    assert code == 2
    done = jobs.Ledger(ledgerPath).done()
    assert 0 < len(done) < 20
    # ... and resuming finishes the rest.
    assert jobs.runTask(Slow(ledgerPath, units = 20), None) == 0
    assert len(jobs.Ledger(ledgerPath).done()) == 20


def test_a_run_over_several_processes_records_each_unit_once(ledgerPath):
    task = Counting(ledgerPath, units = 12)
    assert jobs.runTask(task, None, jobs = 3) == 0
    units = [record['unit'] for record in jobs.Ledger(ledgerPath).records()]
    assert sorted(units, key = int) == [str(index) for index in range(12)]


def test_a_failing_unit_does_not_land_in_the_ledger(ledgerPath):
    """A unit that raised is not finished, so a resumed run must do it again."""
    with pytest.raises(RuntimeError):
        jobs.runTask(Counting(ledgerPath, units = 4, boom = "2"), None)
    assert "2" not in jobs.Ledger(ledgerPath).done()


def test_the_summary_reads_the_ledger_back(ledgerPath, capsys):
    jobs.runTask(Counting(ledgerPath, units = 4), None)
    capsys.readouterr()
    Counting(ledgerPath).summarise(jobs.Ledger(ledgerPath).records(), None)
    assert "[0, 1, 4, 9]" in capsys.readouterr().out
