"""The overnight driver's job table, and the one thing in it that fails silently.

Nothing here runs a job. These are the checks that catch a night wasted before
it is left running: a job naming a task that does not exist, a job with no way
to read its answer back in the morning, and a keep-awake that holds a process
open all night without keeping anything awake.
"""

import argparse
import contextlib
import io

import pytest

import batch
import overnight


def test_every_job_names_a_real_front_door():
    for name, command in overnight.jobs(9).items():
        assert command[0].endswith(".py"), (name, command)
        if command[0] != "batch.py":
            continue
        task = command[1]
        assert task in batch.TASKS, (name, task)


def test_every_batch_job_parses_against_its_task():
    """A job's flags are only checked when it runs, which is after everyone left."""
    for name, command in overnight.jobs(9).items():
        if command[0] != "batch.py":
            continue
        task = batch.TASKS[command[1]]
        parser = argparse.ArgumentParser()
        task.addArguments(parser)
        parser.add_argument("--jobs", type = int, default = 1)
        parser.add_argument("--budget-hours", type = float, dest = "budgetHours")
        args = parser.parse_args(command[2:])
        assert task.ledgerPath(args), name


def test_every_job_says_how_to_read_it_back():
    for name in overnight.jobs(9):
        assert name in overnight.MORNING, name
        summary, aimedAt = overnight.MORNING[name]
        assert summary.startswith("python "), name
        assert aimedAt, name


def test_the_length_15_night_is_three_jobs_that_fit_one_machine():
    """The point of the set: it is asked for by name because it wants the lot."""
    table = overnight.jobs(9)
    night = ['sample15', 'cores15', 'cores13']
    workers = 0
    for name in night:
        command = table[name]
        assert name in table, name
        workers += int(command[command.index("--jobs") + 1])
    assert workers == 15


def test_the_keep_awake_does_not_use_the_spelling_that_silently_fails():
    """`0x80000000 -bor 1` reads as a negative Int32 in PowerShell and throws.

    The first version of this used it, with the result piped to `Out-Null`, so
    the holder sat there for the whole night having kept nothing awake. The
    flags go in as a decimal above Int32's range, which PowerShell reads as
    Int64 and casts exactly.
    """
    snippet = overnight._WINDOWS_STAY_AWAKE
    assert "0x80000000" not in snippet
    assert "[uint32]2147483649" in snippet
    # ES_CONTINUOUS | ES_SYSTEM_REQUIRED, spelled out here so that a change to
    # the literal has to be a deliberate one.
    assert 2147483649 == 0x80000000 | 0x00000001
    assert "exit 1" in snippet, "a failed call must not hold a process open"


@pytest.mark.parametrize("platform", ["win32", "darwin", "linux"])
def test_is_wsl_only_claims_linux(platform, monkeypatch):
    monkeypatch.setattr(overnight.sys, "platform", platform)
    if platform != "linux":
        assert overnight.isWSL() is False


# -- running a night that nobody wrote a job for ---------------------------
#
# The named jobs are one set of parameters each, and the thing a run of nights
# needs is the same tasks at a dozen lengths and widths.  `--run` takes a
# command straight through, so accumulating data does not mean editing this
# file every evening -- and it still gets the keep-awake, the restart and the
# budget, which is the whole reason to go through here at all.

def test_an_ad_hoc_job_is_named_after_the_words_before_its_first_flag():
    assert overnight._adhocName("batch.py cores 16 --max-word 3", {}) == "batch-cores-16"
    assert overnight._adhocName("batch.py sample 14 --count 9000", {}) == "batch-sample-14"
    assert overnight._adhocName("classify.py 10 --resume", {}) == "classify-10"


def test_two_ad_hoc_jobs_with_the_same_name_do_not_share_a_log():
    taken = {"batch-cores-16": []}
    assert overnight._adhocName("batch.py cores 16 --jobs 4",
                                taken) == "batch-cores-16-2"


def test_an_ad_hoc_job_gets_the_night_s_budget_appended():
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        overnight.main(["--dry-run", "--hours", "6",
                        "--run", "batch.py cores 16 --max-word 3 --jobs 8"])
    printed = out.getvalue()
    assert "batch-cores-16" in printed
    assert "--budget-hours 6.0" in printed
    assert "batch.py cores 16 --max-word 3 --jobs 8" in printed


def test_asking_only_for_an_ad_hoc_job_does_not_also_start_the_default_set():
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        overnight.main(["--dry-run", "--run", "batch.py cores 16 --jobs 8"])
    lines = [line for line in out.getvalue().splitlines() if line.startswith("[")]
    assert len(lines) == 1, lines


def test_an_ad_hoc_command_is_split_as_a_shell_would_split_it():
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        overnight.main(["--dry-run", "--run", "batch.py cores 15 --cores 45,504"])
    assert "--cores 45,504" in out.getvalue()
