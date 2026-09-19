#!/usr/bin/env python
"""Long runs, named and resumable, so a machine can be left working on them.

    python batch.py --list                        what there is to run
    python batch.py sample 14 --count 2000        probe 2000 LNAs of length 14
    python batch.py sample 14 --count 2000 --jobs 7 --budget-hours 9
    python batch.py sample 14 --summary           what the ledger says, no work

Every task writes an append-only ledger under `logs/`, one line per finished
unit, and every run resumes from it: the same command again does what is left
and nothing that is done.  `--budget-hours` stops cleanly between units and
exits 2, which is what `overnight.py` restarts on.

The tasks here are the ones that are *only* worth running long.  The three
existing long jobs keep their own front doors, because each has parameters that
do not fit a common shape, and `--list` says what they are.
"""

import argparse
import sys

from quivermutation import jobs
from quivermutation import sampling


class SampleTask(jobs.Task):
    """Draw LNAs of an intractable length and see what the cheap pipeline says.

    Above about n = 13 there are too many LNAs to classify them all -- 208012 at
    n = 13 and 1767263190 at n = 20 -- and the short lengths that *can* be done
    completely are unrepresentative: in a quiver of length 8 every vertex is
    within three arrows of an end.  So the question a long length can be asked
    is a statistical one, and this asks it: what fraction of LNAs of length n
    does the quipu theorem name, what fraction do the moves carry to one it
    names, and what is left over.

    The leftovers are the point.  At n = 9, 10 and 11 they are a handful of
    orbits whose merging is research H-013, and nobody has seen whether their
    rate falls, holds or rises with the length, or whether at n = 16 they take
    shapes that n = 9 has no room for.  Each drawn leftover is recorded with its
    overlap profile, so the shapes can be counted afterwards.

    With `--depth` above 0 each leftover also gets a deduplicated mutation
    search, and what it reaches is recorded.  That is much the most expensive
    part, which is why it is off unless asked for.
    """

    name = 'sample'
    help = "draw LNAs of a long length and probe them"

    def addArguments(self, parser):
        parser.add_argument("length", type = int,
                            help = "the line length to draw from")
        parser.add_argument("--count", type = int, default = 1000,
                            help = "how many LNAs to draw (default 1000)")
        parser.add_argument("--seed", type = int, default = 0,
                            help = "the run's seed; a different seed is a "
                                   "different sample and a different ledger")
        parser.add_argument("--depth", type = int, default = 0,
                            help = "also search this deep from every leftover "
                                   "(default 0, meaning do not search)")
        parser.add_argument("--orbit-limit", type = int, default = 20000,
                            dest = "orbitLimit",
                            help = "how far to walk a move orbit before calling "
                                   "the row a leftover (default 20000)")

    def ledgerPath(self, args):
        # The depth is in the name because a run with a search and a run without
        # are different work on the same draws, and a ledger must not claim a
        # unit was done to a depth it was not.
        return "logs/sample-n{0}-s{1}-d{2}.jsonl".format(
            args.length, args.seed, args.depth)

    def units(self, args):
        return ["{0}/{1}/{2}".format(args.length, args.seed, index)
                for index in range(args.count)]

    def run(self, unit, args):
        index = int(unit.rsplit("/", 1)[1])
        relLengths = sampling.drawFor(args.length, args.seed, index)
        record = sampling.probe(args.length, relLengths,
                                orbitLimit = args.orbitLimit)
        record['index'] = index
        if args.depth > 0 and record['settledBy'] == 'leftover':
            record['search'] = _searchFrom(args.length, relLengths, args.depth)
        return record

    def summarise(self, records, args, out = sys.stdout):
        import collections
        results = [record['result'] for record in records]
        if not results:
            print("nothing in the ledger yet", file = out)
            return
        tally = collections.Counter(result['settledBy'] for result in results)
        total = len(results)
        print("n = {0}, seed {1}: {2} of {3} drawn".format(
            args.length, args.seed, total, args.count), file = out)
        print("  out of {0} LNAs of this length\n".format(
            sampling.countLNAs(args.length)), file = out)
        for settledBy in ('theorem', 'moves', 'leftover'):
            count = tally.get(settledBy, 0)
            print("  {0:<9} {1:6d}  {2:5.1f}%  +- {3:.1f}".format(
                settledBy, count, 100.0 * count / total,
                100.0 * _standardError(count, total)), file = out)

        leftovers = [result for result in results if result['settledBy'] == 'leftover']
        if not leftovers:
            print("\n  no leftovers drawn, so nothing to say about their shape",
                  file = out)
            return
        print("\n  leftovers by maximum relation overlap:", file = out)
        shapes = collections.Counter(result['maxOverlap'] for result in leftovers)
        for overlap in sorted(shapes):
            print("    overlap {0}: {1}".format(overlap, shapes[overlap]), file = out)
        print("\n  leftovers by number of relations:", file = out)
        counts = collections.Counter(result['relations'] for result in leftovers)
        for relations in sorted(counts):
            print("    {0} relations: {1}".format(relations, counts[relations]),
                  file = out)
        print("\n  the first few, by name:", file = out)
        for result in leftovers[:12]:
            print("    {0}  overlap {1}, orbit {2}".format(
                result['name'], result['maxOverlap'], result['orbit']), file = out)


def _searchFrom(length, relLengths, depth):
    """A deduplicated mutation search from one LNA, and what it reached.

    Deduplicated because without it this is unaffordable: the walk revisits the
    same algebra along many routes, 11.4x at n = 9 and depth 6, and the factor
    roughly doubles per level (research E-042).
    """
    from quivermutation import fingerprint, lnaMoves as lm, nakayama as nk, search

    algebra = nk.LinearNakayamaAlgebra(length, list(relLengths))
    visited = fingerprint.Visited()
    collected, hereditary = [], []
    search.mutationSearchDepthFirst(algebra, depth, [], 'sample',
                                    printOutput = False, collected = collected,
                                    collectedHereditary = hereditary,
                                    visited = visited)
    reached = set()
    for pathAlg, _path, _numbering in collected:
        row = lm.asRelLengths(lm._copy(pathAlg), length)
        if row is not None:
            reached.add(''.join(str(value) for value in row))
    return {
        'depth': depth,
        'reached': sorted(reached),
        'hereditary': sorted({quipu for _form, quipu, _path in hereditary if quipu}),
        'walk': visited.summarise(),
    }


def _standardError(count, total):
    """The standard error of a proportion, for reading a sample honestly.

    A run that draws 2000 of 1767263190 LNAs and reports "3.1% leftovers" is
    reporting 3.1% plus or minus 0.4, and the difference between that and a
    later run's 2.6% is not a trend.  Printing the error is cheaper than
    explaining afterwards that it was never there.
    """
    if total <= 0:
        return 0.0
    proportion = count / total
    return (proportion * (1 - proportion) / total) ** 0.5


TASKS = {task.name: task for task in [SampleTask()]}


#: The long jobs that keep their own front door, with the command that runs
#: them.  They are listed here so that `--list` is the whole inventory and a
#: session does not have to remember which script does what.
ELSEWHERE = [
    ("classify", "python classify.py 10 --resume --budget-hours 9",
     "classify a whole length up to derived equivalence"),
    ("merges", "python merges.py 10 --depths 5 6 7 8 --jobs 7 --budget-hours 9",
     "search for mutation paths between the orbits the moves leave over"),
    ("overlaps", "python overlaps.py 10 --free --doubles --no-rules",
     "which relation-overlap configurations a length leaves unplaced"),
    ("discover", "python discover.py --max-arrows 7 --max-width 8",
     "search for new rewrite rules"),
]


def main(argv = None):
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--list", action = "store_true", dest = "listTasks",
                        help = "print the tasks and exit")
    subparsers = parser.add_subparsers(dest = "task")
    for task in TASKS.values():
        sub = subparsers.add_parser(task.name, help = task.help,
                                    description = task.__doc__,
                                    formatter_class = argparse.RawDescriptionHelpFormatter)
        task.addArguments(sub)
        sub.add_argument("--jobs", type = int, default = 1,
                         help = "worker processes (default 1)")
        sub.add_argument("--budget-hours", type = float, default = None,
                         dest = "budgetHours",
                         help = "stop cleanly once this long has passed, "
                                "exiting 2 so a wrapper can restart it")
        sub.add_argument("--summary", action = "store_true",
                         help = "print what the ledger establishes and stop")

    args = parser.parse_args(argv)
    if args.listTasks or args.task is None:
        _printList()
        return 0

    task = TASKS[args.task]
    if args.summary:
        ledger = jobs.Ledger(task.ledgerPath(args))
        task.summarise(ledger.records(), args)
        return 0
    return jobs.runTask(task, args, budgetHours = args.budgetHours, jobs = args.jobs)


def _printList():
    print("Tasks that run here:\n")
    for task in TASKS.values():
        print("  {0:<10} {1}".format(task.name, task.help))
        print("             python batch.py {0} --help".format(task.name))
    print("\nLong jobs with their own front door:\n")
    for name, command, description in ELSEWHERE:
        print("  {0:<10} {1}".format(name, description))
        print("             {0}".format(command))
    print("\n`overnight.py` runs a set of these for a fixed number of hours and "
          "restarts\nwhatever dies.  Every one of them resumes from its own "
          "checkpoint.")


if __name__ == "__main__":
    sys.exit(main())
