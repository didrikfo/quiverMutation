#!/usr/bin/env python
"""Look for local rewrites of an LNA that hold wherever they apply.

    python discover.py                      3 mutations, patterns of <= 5 arrows
    python discover.py --max-steps 4        the same patterns, four mutations
    python discover.py --max-arrows 5 --max-width 7 --jobs 4
    python discover.py --anchor both        rules that need an end of the quiver
    python discover.py --extend             widen the rules already known to tolerate
                                            one relation they do not touch

A *move* is a rewrite of a window of arrows: given the relations inside the
window, replace them with others, by a fixed sequence of mutations at fixed
offsets.  A move is worth having only if it holds at every position of every
length where its window fits, which is what `--verify-lengths` checks.

Discovery plants each candidate pattern of relations in the **middle** of a long
quiver and mutates only near it (research H-007: on a quiver of length 6 or 7
every vertex is within a step of an end, so a rule about the interior cannot be
told apart from one about a boundary).  Two embeddings at different lengths and
offsets are used, and a rewrite is only reported if it recurs across them.

With `--extend` nothing is planted and nothing is mutated speculatively: each
rule already in the table is widened to admit one *spectator* -- a relation
inside its window that the rewrite leaves alone -- and the widened rules are
verified.  That is the cheapest lead there is, because what stops a known rule
from firing is usually a bystander rather than a window that is too small
(research F-023, H-011).  Only the candidates that would fire on an LNA no
search has placed yet are checked.

With `--anchor`, the pattern is planted flush against an end of the quiver
instead, and what comes out is an *anchored* rule -- one stated as holding at
that end and checked only there.  Those are not a curiosity: an isolated pair of
relations overlapping in two or more arrows cannot be pulled apart anywhere in
the interior (research F-021), and the end of the quiver is where it can
(`lnaMoves.endPairCollapseRules`).

The output is the surviving rules as the literal tuples `lnaMoves.VERIFIED_MOVES`
holds, ready to paste in, each with the number of confirmations behind it.
"""

import argparse
import multiprocessing
import sys
import time

from quivermutation import lnaMoves as lm
from quivermutation import overlap as ov


def describeOne(job):
    """One (pattern, length, offset): the rewrites it reaches.

    Top level rather than a closure so multiprocessing can pickle it.
    """
    pattern, length, offset, maxSteps, margin = job
    return lm.discoverLocalMoves([pattern], maxSteps = maxSteps, margin = margin,
                                 embeddings = ((length, offset),), minOccurrences = 1)


def describeOneAnchored(job):
    """One (pattern, anchor, length): the rewrites it reaches against that end."""
    pattern, anchor, length, maxSteps, margin = job
    return lm.discoverAnchoredMoves([pattern], anchor, maxSteps = maxSteps,
                                    margin = margin, lengths = (length,),
                                    minOccurrences = 1)


def spectatorCandidates(rules, maxArrows, maxExtra, lengths):
    """Widened rules that would fire on an LNA the theorem and the table miss.

    Generating them is free; verifying them is not, and most widenings of most
    rules are false.  So the useful ones are picked out first by asking which
    match something still unplaced -- a question about relation lengths only,
    with no mutation in it.
    """
    known = set(rules)
    candidates = set()
    for rule in rules:
        for leftExtra in range(0, maxExtra + 1):
            for rightExtra in range(0, maxExtra + 1 - leftExtra):
                candidates.update(lm.spectatorExtensions(
                    rule, maxArrows = maxArrows,
                    leftExtra = leftExtra, rightExtra = rightExtra))
    candidates -= known

    unplaced = {length: ov.coverage(length, rules)['uncovered'] for length in lengths}
    useful = []
    for candidate in sorted(candidates):
        for length in lengths:
            if any(lm.matchesAt(length, list(relLengths), candidate, windowStart)
                   for relLengths in unplaced[length]
                   for windowStart in lm.windowStartsFor(length, candidate)):
                useful.append(candidate)
                break
    return len(candidates), useful


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
    parser.add_argument("--verify-span", type = int, default = 4, dest = "verifySpan",
                        help = "how many lengths to verify a survivor over, counted from "
                               "the shortest quiver its window fits in (default 4). "
                               "The lengths have to follow the window: a window of 9 "
                               "arrows fits in A_10 at one position only, so a fixed "
                               "range like 7..10 gives such a rule a single confirmation "
                               "from a single length -- and 15 of the 30 window-9 rules "
                               "the first three-mutation run reported that way turned out "
                               "to be false at length 11. E-011.")
    parser.add_argument("--verify-cap", type = int, default = 12, dest = "verifyCap",
                        help = "the longest quiver to verify at (default 12, 58786 LNAs; "
                               "13 is 208012 and 14 is 742900)")
    parser.add_argument("--jobs", type = int, default = 1,
                        help = "parallel worker processes (default 1)")
    parser.add_argument("--known", action = "store_true",
                        help = "also report rewrites already in the table")
    parser.add_argument("--anchor", choices = ("none", "left", "right", "both"),
                        default = "none",
                        help = "plant each pattern against an end of the quiver and "
                               "report anchored rules, rather than in the middle "
                               "(default none). The lengths used are --anchor-lengths.")
    parser.add_argument("--anchor-lengths", default = "11,12", dest = "anchorLengths",
                        help = "quiver lengths to plant against an end at "
                               "(default 11,12)")
    parser.add_argument("--extend", action = "store_true",
                        help = "widen the rules already in the table to tolerate one "
                               "relation they do not touch, instead of searching")
    parser.add_argument("--extend-arrows", type = int, default = 5, dest = "extendArrows",
                        help = "arrows the spectator may have (default 5)")
    parser.add_argument("--extend-margin", type = int, default = 3, dest = "extendMargin",
                        help = "arrows a rule's window may grow by to make room "
                               "(default 3)")
    parser.add_argument("--extend-lengths", default = "7,8,9", dest = "extendLengths",
                        help = "the lengths whose unplaced LNAs a candidate has to "
                               "fire on to be worth verifying (default 7,8,9)")
    parser.add_argument("--max-window", type = int, default = 9, dest = "maxWindow",
                        help = "the widest window to verify at all, since the lengths "
                               "a rule is checked at follow its window (default 9)")
    args = parser.parse_args(argv)

    embeddings = [tuple(int(part) for part in pair.split(":"))
                  for pair in args.embeddings.split(",")]

    if args.extend:
        started = time.time()
        generated, fresh = spectatorCandidates(
            lm.ALL_MOVES, args.extendArrows, args.extendMargin,
            [int(part) for part in args.extendLengths.split(",")])
        fresh = [d for d in fresh if d[0] <= args.maxWindow]
        print("{0} widenings of {1} rules, {2} of them fire on an LNA no search has "
              "placed and fit --max-window, in {3:.0f}s".format(
                  generated, len(lm.ALL_MOVES), len(fresh), time.time() - started))
        pool = multiprocessing.Pool(args.jobs) if args.jobs > 1 else None
        return verifyAndReport(pool, fresh, args)

    patterns = lm.smallPatterns(args.maxRelations, args.maxArrows, args.maxWidth)
    anchors = {"none": (), "left": ("left",), "right": ("right",),
               "both": ("left", "right")}[args.anchor]
    anchorLengths = [int(part) for part in args.anchorLengths.split(",")]
    if anchors:
        jobs = [(pattern, anchor, length, args.maxSteps, args.margin)
                for pattern in patterns for anchor in anchors for length in anchorLengths]
        describe = describeOneAnchored
        print("{0} patterns x {1} ends x {2} lengths = {3} searches, {4} mutations "
              "each".format(len(patterns), len(anchors), len(anchorLengths), len(jobs),
                            args.maxSteps))
    else:
        jobs = [(pattern, length, offset, args.maxSteps, args.margin)
                for pattern in patterns for length, offset in embeddings]
        describe = describeOne
        print("{0} patterns x {1} embeddings = {2} searches, {3} mutations each".format(
            len(patterns), len(embeddings), len(jobs), args.maxSteps))

    pool = multiprocessing.Pool(args.jobs) if args.jobs > 1 else None
    started = time.time()
    seen = {}
    for result in run(pool, describe, jobs):
        for description, places in result.items():
            seen.setdefault(description, []).extend(places)
    print("{0} rewrites described, in {1:.0f}s".format(len(seen), time.time() - started))

    # Recurring across embeddings is the point of planting the pattern twice: a
    # rewrite seen at one embedding only may be an accident of that quiver's ends.
    distinctLengths = len(anchorLengths) if anchors else len(embeddings)
    recurring = {d: places for d, places in seen.items()
                 if len({place[0] for place in places}) >= min(2, distinctLengths)}
    known = set(lm.ALL_MOVES)
    fresh = [d for d in sorted(recurring) if args.known or d not in known]
    print("{0} recur across embeddings, {1} of them new".format(len(recurring), len(fresh)))
    if not fresh:
        return 0

    return verifyAndReport(pool, fresh, args)


def verifyAndReport(pool, fresh, args):
    """Verify each candidate where its window fits, and print the survivors."""
    if not fresh:
        return 0

    def lengthsFor(width):
        """Lengths a window of this width fits in, capped so a run terminates."""
        return [length for length in range(width + 1, width + 1 + args.verifySpan)
                if length <= args.verifyCap]

    started = time.time()
    results = run(pool, verifyOne, [(d, lengthsFor(d[0])) for d in fresh])
    if pool is not None:
        pool.close()

    # Two lengths at least: a rule confirmed at one length only has not been
    # separated from an accident of that quiver's ends, which is the whole point
    # of planting the pattern in the interior.
    survivors = []
    thin = 0
    for description, confirmed, failures in results:
        if failures or confirmed < 2:
            continue
        lengths = lengthsFor(description[0])
        if len(lengths) < 2:
            thin += 1
            continue
        survivors.append((confirmed, description))
    print("{0} verified, in {1:.0f}s{2}".format(
        len(survivors), time.time() - started,
        "" if not thin else "; {0} unverifiable within --verify-cap".format(thin)))
    print()
    for confirmed, description in sorted(survivors, key = lambda pair: -pair[0]):
        print("    {0!r},   # {1} confirmed: {2}".format(
            description, confirmed, lm.formatMove(description)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
