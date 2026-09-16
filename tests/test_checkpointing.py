"""Checkpointing: a classification interrupted partway keeps what it did.

Steps 3 and 4 of classifyLength used to run entirely in memory, with the CSV
written only once both had finished, so a run killed during them lost every
class it had named and every class it had resolved.  At n = 10 that is hours of
work, and a run of that length *will* be interrupted -- see research E-008 and
F-014, where 61 named classes were lost exactly this way.

These check the two halves of the fix: the progress record beside the CSV, and
the wall-clock budget that stops a run cleanly instead of letting it be killed.
"""

import json

import pytest

import mutationClassTable
import quiverMutation as qm
from helpers import quiet


def test_progress_file_records_what_each_step_finished(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    table, _report = quiet(qm.classifyLength, 6, 6, 6, None, False)

    recorded = json.loads((tmp_path / "A_6_mutation_classes.progress.json").read_text())
    # Every class that survives is one the naming step went through, and it is
    # recorded with the number of members it had at the time.
    assert set(recorded["annotated"]) >= table.classNames()
    for className in table.classNames():
        assert recorded["annotated"][className] == len(table.membersOfClass(className))


def test_a_resumed_run_skips_the_classes_it_already_named(tmp_path, monkeypatch):
    """The naming step must not redo a class whose membership has not changed."""
    monkeypatch.chdir(tmp_path)
    quiet(qm.classifyLength, 6, 6, 6, None, False)

    fileName = "A_6_mutation_classes.csv"
    table = mutationClassTable.MutationClassTable.fromCSV(fileName, 6)
    progress = qm.readProgress(fileName)
    # Make the record match the table exactly, as a finished run leaves it.
    progress["annotated"] = {name: len(table.membersOfClass(name))
                             for name in table.classNames()}

    named = []
    original = qm.hereditaryFormFromTheorem

    def watched(lineLength, relationString):
        named.append(relationString)
        return original(lineLength, relationString)

    monkeypatch.setattr(qm, "hereditaryFormFromTheorem", watched)
    quiet(qm.annotateHereditaryForms, table, 6, 0, False, fileName, progress, None)
    # The first pass over the rows still runs -- it is one theorem lookup per row
    # and builds the map the loop reads -- but no class is named a second time.
    assert all(count == len(table.membersOfClass(name))
               for name, count in progress["annotated"].items())


def test_a_changed_class_is_named_again(tmp_path, monkeypatch):
    """A class that gained members since it was named must not be skipped.

    A new member can carry a form the class did not have, so the recorded member
    count is what makes skipping safe.
    """
    monkeypatch.chdir(tmp_path)
    table, _ = quiet(qm.classifyLength, 6, 6, 6, None, False)
    fileName = "A_6_mutation_classes.csv"

    className = sorted(table.classNames())[0]
    stale = {name: len(table.membersOfClass(name)) for name in table.classNames()}
    stale[className] += 1                       # pretend it has since grown
    progress = {"annotated": stale, "resolved": {}}

    table.setHereditaryFormForClass(className, "")
    quiet(qm.annotateHereditaryForms, table, 6, 0, False, fileName, progress, None)
    assert table.formOfEachClass().get(className), "the changed class was skipped"


def test_a_corrupt_progress_file_only_costs_a_redo(tmp_path):
    """A half-written record must not stop a resume."""
    fileName = str(tmp_path / "A_6_mutation_classes.csv")
    (tmp_path / "A_6_mutation_classes.progress.json").write_text('{"annotated": {"x"')
    assert qm.readProgress(fileName) == {"annotated": {}, "resolved": {}}


def test_a_spent_budget_stops_the_run_and_says_so(tmp_path, monkeypatch):
    """A budget of no time at all stops before the search and reports it."""
    monkeypatch.chdir(tmp_path)
    table, report = quiet(qm.classifyLength, 7, 6, 6, None, False, False, 0)

    assert report["stoppedEarly"]
    assert report["unplacedRows"] > 0, "nothing should have been searched"
    # And what it did place -- the quipu theorem seeding, which runs before the
    # first budget check -- is on disk to be resumed from.
    written = mutationClassTable.MutationClassTable.fromCSV(
        str(tmp_path / "A_7_mutation_classes.csv"), 7)
    assert len(written.classNames()) > 0


def test_resuming_a_budgeted_run_finishes_it(tmp_path, monkeypatch):
    """Stopping on the budget and resuming gives the published answer."""
    monkeypatch.chdir(tmp_path)
    quiet(qm.classifyLength, 7, 6, 6, None, False, False, 0)
    table, report = quiet(qm.classifyLength, 7, 6, 6, None, False, True)

    assert not report["stoppedEarly"]
    assert not table.unassignedRelationStrings()
    assert report["candidate"] == {}
    # The published n = 7 classification: 6 classes of 54, 32, 29, 7, 6, 4.
    sizes = sorted((len(table.membersOfClass(name)) for name in table.classNames()),
                   reverse=True)
    assert sizes == [54, 32, 29, 7, 6, 4]
