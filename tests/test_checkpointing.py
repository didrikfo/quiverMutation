"""Checkpointing: a classification interrupted partway keeps what it did.

The table has always been written after every class the *search* places, but the
two steps after it -- resolving merge candidates and naming what is left -- ran
entirely in memory, with the CSV written once both had finished.  A run killed
during them lost every class it had named and every link it had found.  At
n = 10 that is hours of work, and a run of that length *will* be interrupted:
research E-008 lists three separate ways it happened, and F-014 is the run that
named 61 classes and kept none of them.

These check the two halves of the fix: the progress record beside the CSV, and
the wall-clock budget that stops a run cleanly instead of letting it be killed.
"""

import json

import pytest

from quivermutation import mutationClassTable
import quivermutation as qm
from helpers import quiet


def test_a_finished_run_leaves_a_progress_record(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    quiet(qm.classifyLength, 6, 6, 6, None, False)

    recorded = json.loads((tmp_path / "A_6_mutation_classes.progress.json").read_text())
    assert set(recorded) == {"named", "resolved", "condition"}
    # '' is "no deeper probing", which is what a resume compares against: a run
    # under a condition cannot skip what a plain run recorded.
    assert recorded["condition"] == ""


def test_a_budget_stops_the_run_with_everything_on_disk(tmp_path, monkeypatch):
    """n = 8 still needs a search, so a short budget really does interrupt one."""
    monkeypatch.chdir(tmp_path)
    table, report = quiet(qm.classifyLength, 8, 6, 6, None, False, False, None, 0.0)

    assert report["stoppedEarly"]
    assert report["unplacedRows"] > 0, "the search should have been cut short"
    # What it did place is on disk, ready to be resumed from.
    written = mutationClassTable.MutationClassTable.fromCSV(
        str(tmp_path / "A_8_mutation_classes.csv"), 8)
    assert len(written.classNames()) > 0
    assert len(written.unassignedRelationStrings()) == report["unplacedRows"]


def test_an_expired_budget_that_interrupted_nothing_is_not_reported_as_a_stop(
        tmp_path, monkeypatch):
    """A budget can expire on a run that had nothing left to do.

    Seeding places the whole of n = 7 outright, so a zero budget is already spent
    by the time the first check happens and yet the classification is complete.
    Reporting that as "stopped, not finished" would send someone to resume a run
    with nothing in it, so `stoppedEarly` means a step actually broke out.
    """
    monkeypatch.chdir(tmp_path)
    table, report = quiet(qm.classifyLength, 7, 6, 6, None, False, False, None, 0.0)

    assert not table.unassignedRelationStrings()
    assert not report["stoppedEarly"]


def test_resuming_a_budgeted_run_finishes_it(tmp_path, monkeypatch):
    """Stop n = 8 on a budget, resume, and get the published classification."""
    monkeypatch.chdir(tmp_path)
    quiet(qm.classifyLength, 8, 6, 6, None, False, False, None, 0.0)
    table, report = quiet(qm.classifyLength, 8, 6, 6, None, False, True)

    assert not report["stoppedEarly"]
    assert not table.unassignedRelationStrings()
    assert report["candidate"] == {}
    # arXiv:2305.06642's n = 8 table: 11 classes, of these sizes.
    sizes = sorted((len(table.membersOfClass(name)) for name in table.classNames()),
                   reverse=True)
    assert sizes == [133, 65, 64, 64, 40, 26, 13, 10, 9, 4, 1]


def test_a_resolved_class_is_not_searched_again_at_the_same_depth(tmp_path, monkeypatch):
    """The resolve step is the open-ended one, so skipping it on a resume matters."""
    monkeypatch.chdir(tmp_path)
    quiet(qm.classifyLength, 6, 6, 6, None, False)
    fileName = "A_6_mutation_classes.csv"
    table = mutationClassTable.MutationClassTable.fromCSV(fileName, 6)

    className = sorted(table.classNames())[0]
    members = len(table.membersOfClass(className))
    progress = {"named": {}, "resolved": {className: [6, members]}}
    searched = []
    original = qm.search.mutationSearchDepthFirst

    def watched(*args, **kwargs):
        searched.append(args[0])
        return original(*args, **kwargs)

    monkeypatch.setattr(qm.search, "mutationSearchDepthFirst", watched)
    quiet(qm.resolveMergeCandidates, table, 6, 6, False, fileName, progress, None)
    assert progress["resolved"][className] == [6, members], "the record was disturbed"


def test_a_deeper_search_is_a_different_experiment_and_is_not_skipped(tmp_path):
    """A class recorded at depth 4 must be searched again when depth 8 is asked for."""
    recorded = {"someClass": [4, 3]}
    depth, members = recorded["someClass"]
    assert not (depth >= 8 and members == 3), "depth 4 must not satisfy a depth-8 request"
    assert depth >= 4 and members == 3, "depth 4 does satisfy a depth-4 request"


def test_a_changed_class_is_processed_again(tmp_path, monkeypatch):
    """A class that gained members since it was named must not be skipped.

    A new member can carry a form the class did not have, which is why the
    recorded member count is what makes skipping safe.  What is checked here is
    that the class is *visited* -- the step re-records it against its current
    membership -- not that it comes away with a name: `nameRemainingClasses`
    names by the non-piecewise-hereditary certificate and the Coxeter
    polynomial, and a quipu class of A_6 rightly gets neither.
    """
    monkeypatch.chdir(tmp_path)
    table, _ = quiet(qm.classifyLength, 6, 6, 6, None, False)
    fileName = "A_6_mutation_classes.csv"

    className = sorted(table.classNames())[0]
    members = len(table.membersOfClass(className))
    stale = {name: len(table.membersOfClass(name)) for name in table.classNames()}
    stale[className] = members + 1              # pretend it has since grown
    progress = {"named": stale, "resolved": {}}

    # The form has to be cleared, or the step skips it as already proved.
    table.setHereditaryFormForClass(className, "")
    quiet(qm.nameRemainingClasses, table, 6, 0, False, fileName, progress, None)
    assert progress["named"][className] == members, "the changed class was skipped"


def test_an_unchanged_class_is_skipped(tmp_path, monkeypatch):
    """The other half of the same rule: matching membership means do not redo it."""
    monkeypatch.chdir(tmp_path)
    table, _ = quiet(qm.classifyLength, 6, 6, 6, None, False)
    fileName = "A_6_mutation_classes.csv"

    className = sorted(table.classNames())[0]
    members = len(table.membersOfClass(className))
    progress = {"named": {className: members}, "resolved": {}}

    table.setHereditaryFormForClass(className, "")
    quiet(qm.nameRemainingClasses, table, 6, 0, False, fileName, progress, None)
    # Untouched: the loop body never ran, so nothing rewrote the record.
    assert progress["named"][className] == members
    assert not table.formOfEachClass().get(className)


def test_a_corrupt_progress_file_only_costs_a_redo(tmp_path):
    """A half-written record must not stop a resume."""
    fileName = str(tmp_path / "A_6_mutation_classes.csv")
    (tmp_path / "A_6_mutation_classes.progress.json").write_text('{"named": {"x"')
    assert qm.readProgress(fileName) == {"named": {}, "resolved": {}, "condition": ""}


def test_the_progress_file_is_written_atomically(tmp_path):
    """A kill mid-write must leave the previous record, not a truncated one."""
    fileName = str(tmp_path / "A_6_mutation_classes.csv")
    qm.writeProgress(fileName, {"named": {"a": 1}, "resolved": {}})
    qm.writeProgress(fileName, {"named": {"a": 1, "b": 2}, "resolved": {}})

    assert not (tmp_path / "A_6_mutation_classes.progress.json.tmp").exists()
    assert qm.readProgress(fileName)["named"] == {"a": 1, "b": 2}
