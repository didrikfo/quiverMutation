#!/usr/bin/env python
"""Look for local rewrites of an LNA that hold wherever they apply.

    python discover.py                      3 mutations, patterns of <= 5 arrows
    python discover.py --max-steps 4        the same patterns, four mutations
    python discover.py --max-arrows 5 --max-width 7 --jobs 4

A *move* is a rewrite of a window of arrows: given the relations inside the
window, replace them with others, by a fixed sequence of mutations at fixed
offsets.  A move is worth having only if it holds at every position of every
length where its window fits, which is what `--verify-lengths` checks.

Discovery plants each candidate pattern of relations in the **middle** of a long
quiver and mutates only near it (research H-007: on a quiver of length 6 or 7
every vertex is within a step of an end, so a rule about the interior cannot be
told apart from one about a boundary).  Two embeddings at different lengths and
offsets are used, and a rewrite is only reported if it recurs across them.

The output is the surviving rules as the literal tuples `lnaMoves.VERIFIED_MOVES`
holds, ready to paste in, each with the number of confirmations behind it.
"""

import argparse
import multiprocessing
import sys
import time

from quivermutation import lnaMoves as lm


def describeOne(job):
    """One (pattern, length, offset): the rewrites it reaches.

    Top level rather than a closure so multiprocessing can pickle it.
    """
    pattern, length, offset, maxSteps, margin = job
    return lm.discoverLocalMoves([pattern], maxSteps = maxSteps, margin = margin,
                                 embeddings = ((length, offset),), minOccurrences = 1)


def verifyOne(job):
    description, lengths = job
    confirmed, failures = lm.verifyMove(description, lengths)
    return description, confirmed, failures


def run(pool, function, jobs):
    """Map, in parallel when a pool was made, keeping the results in order."""
    if pool is None:
        return [function(job) for job in jobs]
    return pool.map(function, jobs, chunksize = 1)


def main(argv = None):
    parser = argparse.ArgumentParser(
        description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--max-steps", type = int, default = 3, dest = "maxSteps",
                        help = "mutations per sequence (default 3). Each step multiplies "
                               "the search; 4 is the open question of H-008.")
    parser.add_argument("--margin", type = int, default = 3,
                        help = "how far from the planted pattern a mutation may be "
                               "(default 3)")
    parser.add_argument("--max-relations", type = int, default = 3, dest = "maxRelations",
                        help = "relations per pattern (default 3)")
    parser.add_argument("--max-arrows", type = int, default = 4, dest = "maxArrows",
                        help = "arrows per relation in a pattern (default 4)")
    parser.add_argument("--max-width", type = int, default = 5, dest = "maxWidth",
                        help = "arrows a pattern may span (default 5)")
    parser.add_argument("--embeddings", default = "13:4,14:5",
                        help = "length:offset pairs to plant each pattern at "
                               "(default 13:4,14:5)")
    parser.add_argument("--verify-lengths", default = "7,8,9,10", dest = "verifyLengths",
                        help = "lengths to verify a survivor over (default 7,8,9,10). "
                               "A wide rule found at one length is thin evidence; E-010 "
                               "is the run that learned that.")
    parser.add_argument("--jobs", type = int, default = 1,
                        help = "parallel worker processes (default 1)")
    parser.add_argument("--known", action = "store_true",
                        help = "also report rewrites already in VERIFIED_MOVES")
    args = parser.parse_args(argv)

    embeddings = [tuple(int(part) for part in pair.split(":"))
                  for pair in args.embeddings.split(",")]
    verifyLengths = [int(part) for part in args.verifyLengths.split(",")]

    patterns = lm.smallPatterns(args.maxRelations, args.maxArrows, args.maxWidth)
    jobs = [(pattern, length, offset, args.maxSteps, args.margin)
            for pattern in patterns for length, offset in embeddings]
    print("{0} patterns x {1} embeddings = {2} searches, {3} mutations each".format(
        len(patterns), len(embeddings), len(jobs), args.maxSteps))

    pool = multiprocessing.Pool(args.jobs) if args.jobs > 1 else None
    started = time.time()
    seen = {}
    for result in run(pool, describeOne, jobs):
        for description, places in result.items():
            seen.setdefault(description, []).extend(places)
    print("{0} rewrites described, in {1:.0f}s".format(len(seen), time.time() - started))

    # Recurring across embeddings is the point of planting the pattern twice: a
    # rewrite seen at one embedding only may be an accident of that quiver's ends.
    recurring = {d: places for d, places in seen.items()
                 if len({place[0] for place in places}) >= min(2, len(embeddings))}
    known = set(lm.VERIFIED_MOVES)
    fresh = [d for d in sorted(recurring) if args.known or d not in known]
    print("{0} recur across embeddings, {1} of them new".format(len(recurring), len(fresh)))
    if not fresh:
        return 0

    started = time.time()
    results = run(pool, verifyOne, [(d, verifyLengths) for d in fresh])
    if pool is not None:
        pool.close()

    survivors = []
    for description, confirmed, failures in results:
        if failures or not confirmed:
            continue
        survivors.append((confirmed, description))
    print("{0} verified over lengths {1}, in {2:.0f}s".format(
        len(survivors), verifyLengths, time.time() - started))
    print()
    for confirmed, description in sorted(survivors, key = lambda pair: -pair[0]):
        print("    {0!r},   # {1} confirmed: {2}".format(
            description, confirmed, lm.formatMove(description)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
