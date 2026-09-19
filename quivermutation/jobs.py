"""Long runs that survive being stopped: a ledger, a unit of work, a driver.

`classify.py` and `merges.py` each grew their own checkpoint, their own
`--resume`, their own `--budget-hours` and their own summary reader, and each
did it slightly differently.  This is that machinery once, so a new long job is
a `Task` and not another 150 lines of file handling.

The shape it imposes is the one those two arrived at independently, because it
is what makes a run leaveable:

* **Work is divided into units, each named by a string.**  A unit is whatever is
  worth not redoing -- one classified class, one searched orbit member, one
  drawn LNA.  The name is the identity: a ledger with a unit's name in it means
  that unit is finished, whoever finished it and whenever.
* **A finished unit is appended to a JSONL ledger immediately.**  Not at the
  end, not every hundred: a run that is killed loses only what was in flight.
  Append-only means a crash mid-write costs the last line and nothing before it,
  and the reader skips a line it cannot parse for exactly that reason.
* **A run is resumable by default and idempotent.**  Starting the same command
  again does the units the ledger does not have and nothing else.
* **A budget stops it cleanly, between units, with everything on disk.**  Exit
  code 2, which is the convention `overnight.py` already restarts on: 0 and 1
  are finished answers, 2 is "out of time, run me again".

What this deliberately does *not* do is decide anything about the mathematics.
A task says what its units are and how to do one; the ledger holds what came
back; the summary is the task's own reading of its records.  `batch.py` is the
command line over the registry.
"""

import json
import multiprocessing
import os
import sys
import time


class Ledger:
    """An append-only JSONL record of finished units, keyed by unit name.

    Every line is one finished unit: `{"unit": ..., "at": ..., "result": {...}}`.
    A line that will not parse is skipped rather than raising -- the only way to
    get one is a process killed mid-append, and the right response to that is to
    redo the unit, not to refuse to start.
    """

    def __init__(self, path):
        self.path = path
        directory = os.path.dirname(path)
        if directory:
            os.makedirs(directory, exist_ok = True)
        self._handle = None

    def records(self):
        """Every finished unit's record, in the order they were appended."""
        if not os.path.exists(self.path):
            return []
        found = []
        with open(self.path) as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    found.append(json.loads(line))
                except ValueError:
                    continue
        return found

    def done(self):
        """The names of the units the ledger says are finished."""
        return {record['unit'] for record in self.records() if 'unit' in record}

    def record(self, unit, result):
        """Append one finished unit, and flush it, so a kill does not lose it."""
        if self._handle is None:
            self._handle = open(self.path, "a")
        self._handle.write(json.dumps({
            'unit': unit,
            'at': time.strftime("%Y-%m-%dT%H:%M:%S"),
            'result': result,
        }, sort_keys = True) + "\n")
        self._handle.flush()
        os.fsync(self._handle.fileno())

    def close(self):
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def __enter__(self):
        return self

    def __exit__(self, *_exception):
        self.close()
        return False


class Task:
    """One kind of long run: what its units are, how to do one, how to read them.

    Subclasses set `name` and `help`, and implement

    * `addArguments(parser)` -- the run's parameters, as argparse options;
    * `ledgerPath(args)` -- where this parameter set's ledger lives.  Two runs
      whose parameters mean different work must not share a ledger, and two runs
      that mean the same work must;
    * `units(args)` -- the unit names, in the order to do them, as a list or a
      generator.  Cheap: it is called before any work is done, and for a run of
      a million draws it should not build a million of anything else;
    * `run(unit, args)` -- do one unit, return a JSON-serialisable dict.  It must
      be a plain function of its arguments, because it is what a worker process
      is handed;
    * `summarise(records, args, out)` -- print what the ledger establishes.

    `run` is called in a worker process when `--jobs` is more than one, so it
    must not rely on anything the parent set up after import.
    """

    name = None
    help = ""

    def addArguments(self, parser):
        pass

    def ledgerPath(self, args):
        raise NotImplementedError

    def units(self, args):
        raise NotImplementedError

    def run(self, unit, args):
        raise NotImplementedError

    def summarise(self, records, args, out = None):
        out = sys.stdout if out is None else out
        print("{0}: {1} units finished".format(self.name, len(records)), file = out)


def _worker(payload):
    """Unpack what a pool can carry and do one unit.  Top level, so it pickles."""
    task, unit, args = payload
    started = time.time()
    result = task.run(unit, args)
    if isinstance(result, dict):
        result.setdefault('seconds', round(time.time() - started, 3))
    return unit, result


def runTask(task, args, budgetHours = None, jobs = 1, out = None):
    """Do the task's outstanding units, recording each as it finishes.

    Returns 0 when every unit is done and 2 when the budget stopped it, which is
    the convention `overnight.py` restarts on.  Units already in the ledger are
    skipped without being built, so a resumed run costs a read of the ledger and
    nothing else for the work it is not redoing.

    With `jobs` above 1 the units are spread over a process pool and recorded as
    they come back, so the ledger's order is completion order rather than the
    order `units` produced -- which is why a unit's name, and not its position,
    is its identity.

    `out` defaults to `sys.stdout` **at call time**, not at import time: a
    default argument is evaluated once when the module loads, so binding the
    stream there would ignore any later redirection -- which is how a caller
    capturing the progress, a test among them, would silently get nothing.
    """
    out = sys.stdout if out is None else out
    ledger = Ledger(task.ledgerPath(args))
    finished = ledger.done()
    outstanding = [unit for unit in task.units(args) if unit not in finished]
    if finished:
        print("resuming: {0} units already done, {1} to go".format(
            len(finished), len(outstanding)), file = out, flush = True)
    if not outstanding:
        print("nothing to do: every unit is in {0}".format(ledger.path),
              file = out, flush = True)
        return 0

    deadline = None if budgetHours is None else time.time() + budgetHours * 3600
    done = 0
    started = time.time()

    def report(unit):
        elapsed = time.time() - started
        rate = done / elapsed if elapsed else 0.0
        left = (len(outstanding) - done) / rate if rate else float('inf')
        print("  [{0}/{1}] {2}  ({3:.1f}/min, about {4} left)".format(
            done, len(outstanding), unit, rate * 60,
            "unknown" if left == float('inf') else _duration(left)),
            file = out, flush = True)

    stopped = False
    with ledger:
        if jobs <= 1:
            for unit in outstanding:
                if deadline is not None and time.time() > deadline:
                    stopped = True
                    break
                _unit, result = _worker((task, unit, args))
                ledger.record(_unit, result)
                done += 1
                report(_unit)
        else:
            payloads = [(task, unit, args) for unit in outstanding]
            with multiprocessing.Pool(jobs) as pool:
                for unit, result in pool.imap_unordered(_worker, payloads):
                    ledger.record(unit, result)
                    done += 1
                    report(unit)
                    if deadline is not None and time.time() > deadline:
                        stopped = True
                        pool.terminate()
                        break

    print("", file = out)
    task.summarise(Ledger(task.ledgerPath(args)).records(), args, out)
    if stopped:
        print("\nstopped on the budget with {0} units left; "
              "rerun the same command to continue".format(len(outstanding) - done),
              file = out, flush = True)
        return 2
    return 0


def _duration(seconds):
    """A rough human reading of a number of seconds."""
    if seconds < 90:
        return "{0:.0f}s".format(seconds)
    if seconds < 90 * 60:
        return "{0:.0f}m".format(seconds / 60)
    return "{0:.1f}h".format(seconds / 3600)
