#!/usr/bin/env python
"""What the cheap steps cannot place, and which relation patterns are to blame.

    python unplaced.py 9
    python unplaced.py 9 --patterns 12

Seeding from the quipu theorem and expanding along the verified move orbits
places part of the table for free; the rest has to be searched, and searching is
what makes a long classification long.  This measures the rest.

The obstruction is always the same shape.  The theorem covers an LNA whose
consecutive relations overlap in at most one arrow ("almost separate"), so every
row it misses contains at least one *overlapping run* -- a maximal group of
relations chained by overlaps of two or more arrows.  This reports the runs that
actually occur, by how many unplaced rows carry them, which is the list rule
discovery should be aimed at: a rewrite for a common run unlocks every row
carrying it, while a rewrite for a pattern that never occurs unlocks nothing.

This is idea 19 of NOTES.md and hypothesis H-003 -- that the rules found so far
mostly keep an LNA inside the set the theorem already covers, and that the wins
are in the heavily overlapping rows.  It needs no mutation search and runs in
seconds.
"""

import argparse
import collections
import sys

import lnaMoves
import mutationClassTable
import quiverMutation as qm


def overlappingRuns(relLengths):
    """The maximal runs of relations chained by overlaps of two or more arrows.

    Each run is returned as a pattern in the shape `lnaMoves` uses -- a tuple of
    (start, arrows) normalised so the leftmost relation starts at 1 -- so it can
    be handed straight to `discoverLocalMoves` as something to plant.

    A single relation is never a run: one relation on its own is almost
    separate, so it is not what blocks a row.
    """
    relations = [(start, arrows)
                 for start, arrows in enumerate(relLengths, start = 1) if arrows]
    runs = []
    current = []
    for earlier, later in zip(relations, relations[1:]):
        overlap = earlier[0] + earlier[1] - later[0]
        if overlap >= 2:
            if not current:
                current = [earlier]
            current.append(later)
        elif current:
            runs.append(current)
            current = []
    if current:
        runs.append(current)
    return [tuple((start - run[0][0] + 1, arrows) for start, arrows in run)
            for run in runs]


def unplacedRows(length, expandByMoves = True):
    """The LNAs that seeding, and optionally the move orbits, leave unplaced."""
    table = mutationClassTable.MutationClassTable.forLength(
        length,
        [qm.relSetToString(relSet)
         for relSet in qm.generateAllPossibleLineRelations(length)],
    )
    qm.seedTableFromQuipuTheorem(table, length, printOutput = False,
                                 expandByMoves = expandByMoves)
    return table, table.unassignedRelationStrings()


def main(argv = None):
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("length", type = int, help = "number of vertices in the line quiver")
    parser.add_argument("--patterns", type = int, default = 15,
                        help = "how many of the commonest runs to list (default 15)")
    parser.add_argument("--no-moves", action = "store_true",
                        help = "seed from the theorem only, without expanding along "
                               "move orbits, to show what the moves are worth")
    args = parser.parse_args(argv)

    length = args.length
    table, unplaced = unplacedRows(length, expandByMoves = not args.no_moves)
    _, theoremOnly = unplacedRows(length, expandByMoves = False)
    total = len(table)
    placed = total - len(unplaced)

    print("A_{0}: {1} LNAs".format(length, total))
    print("  placed by the quipu theorem alone   {0:>6}  ({1:.0f}%)".format(
        total - len(theoremOnly), 100 * (total - len(theoremOnly)) / total))
    if not args.no_moves:
        print("  placed by theorem + move orbits     {0:>6}  ({1:.0f}%)".format(
            placed, 100 * placed / total))
    print("  still needing a search              {0:>6}  ({1:.0f}%)".format(
        len(unplaced), 100 * len(unplaced) / total))
    if not unplaced:
        return 0

    runCounts = collections.Counter()
    rowsWithRun = collections.Counter()
    runsPerRow = collections.Counter()
    for relationString in unplaced:
        relLengths = qm.relationStringToLineRelLengths(length, relationString)
        runs = overlappingRuns(relLengths)
        runsPerRow[len(runs)] += 1
        for run in runs:
            runCounts[run] += 1
        for run in set(runs):
            rowsWithRun[run] += 1

    print()
    print("Overlapping runs per unplaced row:")
    for count, rows in sorted(runsPerRow.items()):
        note = "  <- not an overlap at all; the moves missed these" if count == 0 else ""
        print("  {0} run(s): {1:>6} rows{2}".format(count, rows, note))

    print()
    print("The commonest overlapping runs, by how many unplaced rows carry one:")
    print("  {0:<34} {1:>7} {2:>9}".format("run (start:arrows)", "rows", "share"))
    for run, rows in rowsWithRun.most_common(args.patterns):
        shape = " ".join("({0}:{1})".format(start, arrows) for start, arrows in run)
        print("  {0:<34} {1:>7} {2:>8.1f}%".format(
            shape, rows, 100 * rows / len(unplaced)))

    covered = set()
    ordered = [run for run, _ in rowsWithRun.most_common()]
    for takeTop in (1, 3, 5, 10, len(ordered)):
        if takeTop > len(ordered):
            continue
        wanted = set(ordered[:takeTop])
        covered = sum(
            1 for relationString in unplaced
            if wanted & set(overlappingRuns(
                qm.relationStringToLineRelLengths(length, relationString))))
        print()
        print("A rule for each of the top {0} run(s) would touch {1} of {2} "
              "unplaced rows ({3:.0f}%).".format(
                  takeTop, covered, len(unplaced), 100 * covered / len(unplaced)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
