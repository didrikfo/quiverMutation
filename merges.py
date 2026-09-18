#!/usr/bin/env python
"""Search for mutation paths between the orbits the moves leave over.

    python merges.py 10                                  depths 5, 6, 7, all cores but one
    python merges.py 10 --depths 5 6 --jobs 7 --budget-hours 9
    python merges.py 11 --depths 4 5 6 --jobs 7
    python merges.py 10 --summary                        what the checkpoint says, no search

Seeding from the quipu theorem, the free move and the double mutation of
arXiv:2310.08346 place every LNA that is in a quipu class at n = 9, 10 and 11,
and leave the rest in a handful of orbits (research F-032). What those orbits
do not say is whether two of them sharing a Coxeter polynomial are one derived
class -- research H-013. This settles that the only way available without a
new invariant: a mutation search out of every member of every such orbit, and
out of its relation dual, looking for a member of a sibling.

The search deepens iteratively -- every member at the first depth, then every
member at the next -- because the cost grows about fivefold per level while
most links, where there are any, are expected shallow. The orbits are closed
under the relation dual first, which is free. A polynomial group
whose orbits have all merged is not searched further.

**A negative answer is a lower bound on depth, not a separation.** Two orbits
that never meet are only known not to meet within the depth searched.

Every finished search is appended to a JSONL checkpoint as it completes, so a
killed run loses only what was in flight and `--resume` is the default. A link
between orbits with *different* Coxeter polynomials is impossible and is
printed as an ALARM: it would mean a bug in the moves or the search.

On Windows the process asks the system not to sleep while it runs
(`SetThreadExecutionState`, the same request a media player makes). Closing a
laptop's lid can suspend it regardless.
"""

import argparse
import collections
import json
import multiprocessing
import os
import sys
import time

from quivermutation import doubleMutation as dm
from quivermutation import freeMoves as fm
from quivermutation import invariants
from quivermutation import lines
from quivermutation import lnaMoves as lm
from quivermutation import nakayama as nk
from quivermutation import overlap as ov
from quivermutation import piecewiseHereditary as ph
from quivermutation import search


def keepAwake():
    if sys.platform != 'win32':
        return
    import ctypes
    ES_CONTINUOUS, ES_SYSTEM_REQUIRED = 0x80000000, 0x00000001
    ctypes.windll.kernel32.SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED)


def coxeter(length, relLengths):
    import sympy
    poly = lm._quiet(invariants.coxeterPoly, nk.LinearNakayamaAlgebra(length, list(relLengths)))
    return str(sympy.factor(poly.as_expr()))


def _coxeterTask(args):
    return coxeter(*args)


def relationDual(length, relLengths):
    return dm.relLengthsOf(length, dm.dualIntervals(length, dm.intervalsOf(relLengths)))


def leftoverGroups(length, pool):
    """The orbits no seed reaches, closed under the relation dual, grouped by polynomial.

    The relation dual of an LNA is derived equivalent to it, and neither the free
    move nor the double mutation joins the two -- a first run at n = 10 spent its
    whole depth-3 pass rediscovering exactly those links, because the search
    starts from each member's dual.
    """
    lnas, byRoot = fm.derivedOrbits(length, rules = [], free = True, doubles = True)
    rootOf = {m: root for root, members in byRoot.items() for m in members}
    joined = Unions()
    for lna in lnas:
        joined.union(rootOf[lna], rootOf[relationDual(length, lna)])
    merged = collections.defaultdict(list)
    for root, members in byRoot.items():
        merged[joined.find(root)].extend(members)
    orbits = merged
    orbitOf = {}
    left = []
    for members in orbits.values():
        name = min(lines.className(m) for m in members)
        for member in members:
            orbitOf[member] = name
        if not any(ov.isAlmostSeparate(length, m) for m in members):
            left.append((name, sorted(members, key = lambda m: (sum(1 for a in m if a), m))))
    polys = pool.map(_coxeterTask, [(length, members[0]) for _, members in left])
    cache = {}
    groups = collections.defaultdict(list)
    for (name, members), poly in zip(left, polys):
        certified = sum(1 for m in members
                        if ph.isNotPiecewiseHereditary(length, list(m))
                        or ph.notPiecewiseHereditaryByDeletion(length, list(m), cache))
        groups[poly].append({'name': name, 'members': members, 'certified': certified})
    return orbitOf, dict(groups)


def searchFrom(args):
    """Every LNA a depth-bounded search reaches from a member and from its dual."""
    length, relLengths, depth = args
    started = time.time()
    relationString = nk.LinearNakayamaAlgebra(length, list(relLengths)).relationString()
    reached = set()
    for startPoint in search.memberAndItsDual(length, relationString):
        collected = []
        search.mutationSearchDepthFirst(startPoint, depth, [], 'merges',
                                        printOutput = False, collected = collected)
        for pathAlg, _path, _numbering in collected:
            row = lm.asRelLengths(lm._copy(pathAlg), length)
            if row is not None:
                reached.add(tuple(row))
    return tuple(relLengths), depth, sorted(reached), time.time() - started


class Unions:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def summarise(length, groups, unions, searched, out = sys.stdout):
    print("\nA_{0}: {1} polynomial groups left over, {2} orbits".format(
        length, len(groups), sum(len(g) for g in groups.values())), file = out)
    total = 0
    for poly, orbits in sorted(groups.items(), key = lambda kv: -sum(len(o['members']) for o in kv[1])):
        classes = collections.defaultdict(list)
        for orbit in orbits:
            classes[unions.find(orbit['name'])].append(orbit)
        total += len(classes)
        print("  {0}".format(poly), file = out)
        print("    {0} orbit(s) -> at most {1} class(es)".format(len(orbits), len(classes)), file = out)
        for orbit in orbits:
            depths = [searched.get(m, 0) for m in orbit['members']]
            print("      {0:>12s}  {1:4d} members, {2:4d} certified not p.h., merged into {3}, "
                  "every member searched to depth {4}".format(
                      orbit['name'], len(orbit['members']), orbit['certified'],
                      unions.find(orbit['name']), min(depths)), file = out)
    print("  at most {0} non-quipu classes (each polynomial group is at least one)".format(total),
          file = out)
    return total


def main(argv = None):
    parser = argparse.ArgumentParser(description = __doc__,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("length", type = int)
    parser.add_argument("--depths", type = int, nargs = "+", default = [5, 6, 7])
    parser.add_argument("--jobs", type = int, default = max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument("--budget-hours", type = float, default = None, dest = "budgetHours")
    parser.add_argument("--checkpoint", default = None,
                        help = "JSONL file (default logs/merges-<n>.jsonl)")
    parser.add_argument("--all-groups", action = "store_true", dest = "allGroups",
                        help = "also search orbits alone in their polynomial group, "
                               "as a check that they reach nothing seeded")
    parser.add_argument("--summary", action = "store_true",
                        help = "print what the checkpoint establishes and stop")
    args = parser.parse_args(argv)

    keepAwake()
    length = args.length
    checkpoint = args.checkpoint or os.path.join("logs", "merges-{0}.jsonl".format(length))
    os.makedirs(os.path.dirname(checkpoint) or ".", exist_ok = True)
    deadline = None if args.budgetHours is None else time.time() + args.budgetHours * 3600

    with multiprocessing.Pool(args.jobs) as pool:
        started = time.time()
        orbitOf, groups = leftoverGroups(length, pool)
        polyOf = {orbit['name']: poly for poly, orbits in groups.items() for orbit in orbits}
        print("A_{0}: {1} leftover orbits in {2} polynomial groups, computed in {3:.0f}s".format(
            length, len(polyOf), len(groups), time.time() - started), flush = True)

        unions = Unions()
        searched = {}
        if os.path.exists(checkpoint):
            with open(checkpoint) as handle:
                for line in handle:
                    try:
                        record = json.loads(line)
                    except ValueError:
                        continue            # a line cut off by a kill
                    member = tuple(record['member'])
                    searched[member] = max(searched.get(member, 0), record['depth'])
                    for other in record['reachedOrbits']:
                        unions.union(orbitOf[member], other)
            print("resumed {0} searches from {1}".format(len(searched), checkpoint), flush = True)

        if args.summary:
            summarise(length, groups, unions, searched)
            return 0

        def settled(poly):
            return len({unions.find(o['name']) for o in groups[poly]}) == 1

        wanted = [poly for poly, orbits in groups.items() if len(orbits) >= 2 or args.allGroups]
        tasks = collections.deque(
            (poly, member, depth)
            for depth in sorted(args.depths)
            for poly in sorted(wanted, key = lambda p: sum(len(o['members']) for o in groups[p]))
            for orbit in groups[poly]
            for member in orbit['members'])

        inFlight = {}
        stoppedOnBudget = False
        with open(checkpoint, "a") as log:
            while tasks or inFlight:
                while tasks and len(inFlight) < args.jobs * 2:
                    if deadline is not None and time.time() > deadline:
                        stoppedOnBudget = True
                        tasks.clear()
                        break
                    poly, member, depth = tasks.popleft()
                    if searched.get(member, 0) >= depth or (settled(poly) and not args.allGroups):
                        continue
                    inFlight[(member, depth)] = pool.apply_async(searchFrom, ((length, member, depth),))
                done = [key for key, result in inFlight.items() if result.ready()]
                if not done:
                    time.sleep(1)
                    continue
                for key in done:
                    member, depth, reached, seconds = inFlight.pop(key).get()
                    own = orbitOf[member]
                    others = sorted({orbitOf[r] for r in reached if r in orbitOf} - {own})
                    # Three kinds of link, and the first version conflated the
                    # last two by asking `polyOf.get(o) != polyOf[own]`: `polyOf`
                    # holds only the leftover orbits, so *every* link to a seeded
                    # orbit came out as an ALARM and was dropped from the union.
                    # That is the most interesting result this run can produce --
                    # an LNA the quipu theorem misses turning out to be in a
                    # theorem class after all -- and it was the one being hidden.
                    covered = [o for o in others if o not in polyOf]
                    alarms = [o for o in others if o in polyOf and polyOf[o] != polyOf[own]]
                    for other in others:
                        if other not in alarms and other not in covered:
                            unions.union(own, other)
                    searched[member] = max(searched.get(member, 0), depth)
                    log.write(json.dumps({'length': length, 'member': list(member),
                                          'class': lines.className(member), 'orbit': own,
                                          'depth': depth, 'lnasReached': len(reached),
                                          'reachedOrbits': [o for o in others
                                                            if o not in alarms and o not in covered],
                                          'coveredOrbits': covered,
                                          'alarms': alarms, 'seconds': round(seconds, 1),
                                          'at': time.strftime('%Y-%m-%d %H:%M:%S')}) + "\n")
                    log.flush()
                    note = ""
                    if others:
                        note = "  LINK to " + ", ".join(others)
                    if covered:
                        note += ("  COVERED: reaches the seeded orbit(s) " + ", ".join(covered)
                                 + " -- this leftover is in a quipu class; verify the path")
                    if alarms:
                        note += ("  ALARM: different Coxeter polynomial -- " + ", ".join(alarms)
                                 + " -- with the Coxeter guard on this should not happen (F-038)")
                    print("{0} depth {1} from {2} (orbit {3}): {4} LNAs, {5:.0f}s{6}".format(
                        time.strftime('%H:%M:%S'), depth, lines.className(member), own,
                        len(reached), seconds, note), flush = True)

        total = summarise(length, groups, unions, searched)
        summaryPath = os.path.join(os.path.dirname(checkpoint) or ".",
                                   "merges-{0}-summary.txt".format(length))
        with open(summaryPath, "w") as out:
            summarise(length, groups, unions, searched, out)
        print("summary written to", summaryPath)
        if stoppedOnBudget:
            print("stopped on the budget; rerun the same command to continue")
            return 2
        return 0


if __name__ == "__main__":
    sys.exit(main())
