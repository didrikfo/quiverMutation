"""The commands in `OVERNIGHT.md` are the interface; check they still run.

That file exists so that a run of nights needs no code changes, which makes it
the one document whose commands are *meant* to be copied unread at midnight. A
flag renamed here and not there costs a whole night to a usage error, and the
error is not seen until morning. So every command in it is parsed against the
task it names, exactly as `batch.py` would parse it.

Nothing here runs any work: `units` is cheap by contract, because `jobs.runTask`
calls it before doing anything.
"""

import argparse
import os
import re
import shlex

import pytest

import batch


DOC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "OVERNIGHT.md")
WORKER_CAP = 14


def _text():
    with open(DOC, encoding = "utf-8") as handle:
        return handle.read()


def _batchCommands():
    """Every `batch.py ...` command the document offers, however it is wrapped."""
    text = _text()
    found = re.findall(r"--run '([^']+)'", text)
    found += re.findall(r"\.venv/bin/python (batch\.py [^\"]*?)(?=\")", text)
    commands = []
    for command in found:
        if not command.startswith("batch.py"):
            continue
        # One of the commands is inside a shell loop over the lengths. Its
        # flags are what this file is checking; the loop variable and the
        # `done` that closes the loop belong to the shell, so stand a length
        # in for the one and drop the other.
        command = command.split(";")[0].strip()
        command = re.sub(r"\\?\$n\b", "15", command)
        commands.append(command)
    return commands


def _parse(command):
    words = shlex.split(command)
    task = batch.TASKS[words[1]]
    parser = argparse.ArgumentParser(prog = words[1])
    task.addArguments(parser)
    parser.add_argument("--jobs", type = int, default = 1)
    parser.add_argument("--summary", action = "store_true")
    parser.add_argument("--plan", action = "store_true")
    parser.add_argument("--budget-hours", type = float, dest = "budgetHours")
    return task, parser.parse_args(words[2:])


def test_the_document_still_offers_commands():
    # A regex that has quietly stopped matching would make every test below
    # pass by finding nothing, which is the failure mode to rule out first.
    assert len(_batchCommands()) >= 20


@pytest.mark.parametrize("command", _batchCommands())
def test_every_command_in_the_document_parses(command):
    task, args = _parse(command)
    assert task.name in batch.TASKS


@pytest.mark.parametrize("command", _batchCommands())
def test_every_command_in_the_document_names_units_that_exist(command):
    task, args = _parse(command)
    units = task.units(args)
    assert units, "a command that asks for nothing is a wasted night"
    assert task.ledgerPath(args).startswith("logs/")


def test_no_night_asks_for_more_workers_than_the_machine_has():
    for block in _text().split("```bash")[1:]:
        block = block.split("```")[0]
        total = sum(int(count) for count in re.findall(r"--jobs (\d+)", block))
        assert total <= WORKER_CAP, block.strip()[:120]


def test_the_flags_the_document_calls_filters_really_are_filters():
    # The document tells the reader that `--cores`, `--core-limit` and `--count`
    # can be varied between nights without splitting a ledger. If that stopped
    # being true, two nights would silently redo each other's work.
    task = batch.CoresTask()
    parser = argparse.ArgumentParser()
    task.addArguments(parser)
    whole = parser.parse_args(["15", "--max-word", "4"])
    for narrowing in (["--core-limit", "40"], ["--cores", "45"]):
        part = parser.parse_args(["15", "--max-word", "4"] + narrowing)
        assert task.ledgerPath(part) == task.ledgerPath(whole)

    sampler = batch.SampleTask()
    parser = argparse.ArgumentParser()
    sampler.addArguments(parser)
    assert (sampler.ledgerPath(parser.parse_args(["15", "--count", "4000"]))
            == sampler.ledgerPath(parser.parse_args(["15", "--count", "50000"])))


def test_the_flags_the_document_calls_ledger_naming_really_are():
    # And the other direction, which is the expensive one to get wrong: a flag
    # that changes what an answer *means* must not share a ledger with one that
    # means something else.
    task = batch.CoresTask()
    parser = argparse.ArgumentParser()
    task.addArguments(parser)
    base = ["15", "--max-word", "4"]
    ledgers = {task.ledgerPath(parser.parse_args(base))}
    for differing in (["--max-word", "3"], ["--max-arrows", "8"],
                      ["--pair-word", "3"], ["--gaps", "1,2"],
                      ["--orbit-limit", "60000"], ["--join-limit", "40000"]):
        path = task.ledgerPath(parser.parse_args(["15"] + differing
                                                 if differing[0] == "--max-word"
                                                 else base + differing))
        assert path not in ledgers, differing
        ledgers.add(path)

    sampler = batch.SampleTask()
    parser = argparse.ArgumentParser()
    sampler.addArguments(parser)
    paths = {sampler.ledgerPath(parser.parse_args(["15"] + flags))
             for flags in ([], ["--seed", "1"], ["--depth", "4"],
                           ["--orbit-limit", "100000"])}
    assert len(paths) == 4
