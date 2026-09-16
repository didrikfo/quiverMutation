#!/usr/bin/env python
"""How far the quipu theorem and the move rules reach, read by relation overlap.

    python overlaps.py 9                 the coverage at one length
    python overlaps.py 6 7 8 9           several, as one table
    python overlaps.py 9 --cores         what the rows still needing a search look like
    python overlaps.py 9 --floating      the floating rules only, leaving the ends out

The quipu theorem names the class of an LNA whose consecutive relations share at
most one arrow.  Everything else has to be searched for -- unless a move rule
carries it to something the theorem does cover.  Overlap is therefore the
coordinate to measure the gap in, and this prints it: how many LNAs sit at each
maximum overlap, how many of them the theorem plus the move orbits place, and
how many are left.

`--cores` prints the heavily overlapping sub-patterns of what is left, commonest
first.  That is the list discovery should be aimed at (research H-003).
"""

import argparse
import sys

from quivermutation import lnaMoves as lm
from quivermutation import overlap as ov


def formatPattern(pattern):
    return " ".join("({0}:{1})".format(start, arrows) for start, arrows in pattern)


def report(length, rules, showCores, coreLimit):
    result = ov.coverage(length, rules)
    total = len(result['lnas'])
    print("\nA_{0}: {1} LNAs, {2} named by the theorem, {3} covered with the move "
          "orbits ({4:.0%}), {5} left, in {6} orbits".format(
              length, total, len(result['seeded']), len(result['covered']),
              len(result['covered']) / total, len(result['uncovered']),
              len(result['orbits'])))
    print("    max overlap    LNAs   covered   left")
    for overlapValue, (count, covered, left) in result['byMaxOverlap'].items():
        print("       {0:^9d} {1:7d} {2:9d} {3:6d}".format(
            overlapValue, count, covered, left))
    if showCores:
        cores = ov.blockingCores(length, result['uncovered'])
        print("    {0} distinct overlapping runs among the {1} left:".format(
            len(cores), len(result['uncovered'])))
        for pattern, count in cores.most_common(coreLimit):
            print("      {0:6d}  {1}   width {2}".format(
                count, formatPattern(pattern), lm.patternWidth(pattern)))
    return result


def main(argv = None):
    parser = argparse.ArgumentParser(
        description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("lengths", type = int, nargs = "+",
                        help = "the quiver lengths to report on")
    parser.add_argument("--cores", action = "store_true",
                        help = "also print the overlapping runs of what is left")
    parser.add_argument("--core-limit", type = int, default = 12, dest = "coreLimit",
                        help = "how many of those to print (default 12)")
    parser.add_argument("--floating", action = "store_true",
                        help = "use only the rules that hold at every position, "
                               "leaving out the ones anchored to an end")
    args = parser.parse_args(argv)

    rules = lm.VERIFIED_MOVES if args.floating else lm.ALL_MOVES
    print("{0} rules: {1} floating{2}".format(
        len(rules), len(lm.VERIFIED_MOVES),
        "" if args.floating else " + {0} anchored".format(len(lm.ANCHORED_MOVES))))
    for length in args.lengths:
        report(length, rules, args.cores, args.coreLimit)
    return 0


if __name__ == "__main__":
    sys.exit(main())
