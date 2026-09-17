#!/usr/bin/env python
"""Ask what a named relation pattern can be turned into, and how deep it takes.

    python probe.py 1:3,2:3                       three mutations, in the interior
    python probe.py 1:3,2:3 --steps 6             deeper
    python probe.py 1:3,2:3 --margin 8            wider
    python probe.py 1:5,2:7 --anchor left         planted against an end instead
    python probe.py 1:3,2:6 1:5,2:6 --steps 4     several patterns in one run

Discovery answers "what rules are there"; this answers "what can happen to
*this* configuration", which is the question a hypothesis about a particular
pattern needs.  A pattern is planted in a quiver -- in the middle by default,
flush against an end with `--anchor` -- and every admissible mutation sequence
within `--margin` vertices of it is enumerated.  What comes back is the LNAs
reached, grouped by how much their relations overlap.

Overlap is the point of the grouping.  The quipu theorem covers exactly the
LNAs whose consecutive relations share at most one arrow (research F-021), so a
sequence that *lowers* the maximum overlap is a sequence that moves a pattern
towards being named outright.  A run printing "nothing lower" is a negative
result and is worth recording as one: that is how F-022 was established, and
research H-010 is the conjecture that no depth will ever change it.

The quiver is chosen to fit the pattern with room on both sides unless
`--length` says otherwise, so that no end is in reach of a pattern meant to be
in the interior (research H-007).
"""

import argparse
import sys
import time

from quivermutation import lines
from quivermutation import lnaMoves as lm
from quivermutation import overlap as ov


def parsePattern(text):
    """`1:3,2:3` as ((1, 3), (2, 3)) -- start vertex and arrow count each."""
    return tuple(tuple(int(part) for part in relation.split(":"))
                 for relation in text.split(","))


def formatPattern(pattern):
    return " ".join("({0}:{1})".format(start, arrows) for start, arrows in pattern)


def placement(pattern, anchor, length, clearance, margin, allowEnds = False):
    """Where to plant the pattern, and in how long a quiver.

    In the interior the pattern needs `clearance` arrows of empty quiver on each
    side, so that a mutation within `margin` of it cannot see an end.  Against an
    end it needs the clearance on one side only.

    The clearance is raised to the margin where the margin is larger, which is
    the whole point of a wide probe: widening the margin without lengthening the
    quiver would let the mutations reach an end, and an interior probe that
    touches an end is measuring the wrong thing (research H-007, F-022).

    `allowEnds` turns that off and takes `--length` at its word, which is a
    different and stronger question: not "what can happen near this pattern"
    but "what can happen to this LNA at all", with every vertex of the quiver
    mutable.  A negative answer there says more than a negative answer here.
    """
    if not allowEnds:
        # A mutation at vertex v touches the arrows v - 1 and v, so the arrows the
        # sequence can reach start at centreLo - margin - 1.  Leaving margin + 2
        # arrows clear is what puts arrow 1 out of reach; margin + 1 leaves vertex
        # 2 mutable, and mutating there rewrites arrow 1.
        clearance = max(clearance, margin + 2)
    width = lm.patternWidth(pattern)
    if anchor == 'left':
        offset, needed = 0, width + clearance + 1
    elif anchor == 'right':
        offset, needed = None, width + clearance + 1
    else:
        offset, needed = clearance, width + 2 * clearance + 1
    length = max(length or 0, needed)
    if anchor == 'right':
        offset = length - 1 - width
    return length, offset


def probe(pattern, anchor, length, clearance, steps, margin, allowEnds = False):
    """Plant the pattern and enumerate what the mutations near it reach."""
    length, offset = placement(pattern, anchor, length, clearance, margin, allowEnds)
    relLengths = lm.embedPattern(length, pattern, offset)
    if relLengths is None:
        return None
    centreLo = offset + min(start for start, _ in pattern)
    centreHi = offset + max(start + arrows - 1 for start, arrows in pattern) - 1
    started = time.time()
    reached = lm.localMutationSequences(length, relLengths, centreLo, centreHi,
                                        steps, margin)
    return {
        'length': length,
        'offset': offset,
        'arrowsInReach': (max(1, centreLo - margin - 1),
                          min(length - 1, centreHi + margin)),
        'relLengths': relLengths,
        'reached': reached,
        'seconds': time.time() - started,
    }


def report(pattern, anchor, result, showLowered):
    relLengths = result['relLengths']
    start = ov.maxOverlap(relLengths)
    byOverlap = {}
    lowered = []
    for name, sequence in result['reached'].items():
        row = [int(character) for character in name]
        overlap = ov.maxOverlap(row)
        byOverlap[overlap] = byOverlap.get(overlap, 0) + 1
        if overlap < start:
            lowered.append((overlap, name, sequence))

    lowArrow, highArrow = result['arrowsInReach']
    touchesEnd = lowArrow <= 1 or highArrow >= result['length'] - 1
    if anchor is not None:
        where = "the {0} end".format(anchor)
    elif touchesEnd:
        where = "A QUIVER WHOSE END IS IN REACH,"
    else:
        where = "the interior"
    print("\n{0}  in {1} of A_{2} at offset {3}".format(
        formatPattern(pattern), where, result['length'], result['offset']))
    print("    {0}, maximum overlap {1}".format(lines.className(relLengths), start))
    print("    mutations can rewrite the arrows {0} to {1} of {2}".format(
        lowArrow, highArrow, result['length'] - 1))
    print("    {0} LNAs reached in {1:.0f}s; by maximum overlap {2}".format(
        len(result['reached']), result['seconds'], dict(sorted(byOverlap.items()))))
    if not lowered:
        print("    NOTHING LOWER -- the overlap does not come down here")
        return False
    lowered.sort()
    print("    LOWERED to {0}: {1} via {2}".format(*lowered[0]))
    for overlap, name, sequence in lowered[1:showLowered]:
        print("      also {0}: {1} via {2}".format(overlap, name, sequence))
    return True


def main(argv = None):
    parser = argparse.ArgumentParser(
        description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("patterns", nargs = "+",
                        help = "patterns as start:arrows pairs, e.g. 1:3,2:3")
    parser.add_argument("--steps", type = int, default = 3,
                        help = "mutations per sequence (default 3). Each step "
                               "multiplies the search")
    parser.add_argument("--margin", type = int, default = 3,
                        help = "how far from the pattern a mutation may be (default 3)")
    parser.add_argument("--clearance", type = int, default = 4,
                        help = "arrows of empty quiver to leave beside the pattern "
                               "(default 4), so that no end is in reach. Raised to "
                               "--margin + 1 when that is larger, since otherwise a "
                               "wide probe would reach an end")
    parser.add_argument("--length", type = int, default = 0,
                        help = "the quiver length to use, if the one the clearance "
                               "implies is not wanted")
    parser.add_argument("--anchor", choices = ("left", "right"), default = None,
                        help = "plant the pattern flush against this end instead of "
                               "in the middle")
    parser.add_argument("--allow-ends", action = "store_true", dest = "allowEnds",
                        help = "do not lengthen the quiver to keep the ends out of "
                               "reach, so that --length and a large --margin together "
                               "ask what can happen to the LNA at all")
    parser.add_argument("--show-lowered", type = int, default = 4, dest = "showLowered",
                        help = "how many overlap-lowering results to print (default 4)")
    args = parser.parse_args(argv)

    anyLowered = False
    for text in args.patterns:
        pattern = parsePattern(text)
        result = probe(pattern, args.anchor, args.length, args.clearance,
                       args.steps, args.margin, args.allowEnds)
        if result is None:
            print("\n{0}: does not fit, or is not an admissible pattern".format(text))
            continue
        anyLowered |= report(pattern, args.anchor, result, args.showLowered)
    return 0 if anyLowered else 1


if __name__ == "__main__":
    sys.exit(main())
