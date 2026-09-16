"""How much the relations of an LNA overlap, and how far the move table reaches.

The quipu theorem of arXiv:2305.06642 names the class of an LNA whose relations
are *almost separate* -- consecutive relations sharing at most one arrow.  That
condition is exactly a bound on the overlaps computed here, so the overlap
profile is the natural coordinate for asking what the theorem misses and what
the move rules of `lnaMoves` add on top of it.

`coverage` is the measurement research H-003 asks for: of the rows a search
still has to place, what relation patterns do they have?  The answer is read
off the overlap profile -- see research F-021.
"""

import collections

from . import lnaMoves
from . import nakayama
from . import quipuForms


def overlapProfile(relLengths):
    """The arrows each consecutive pair of relations shares, left to right.

    Relation (start, arrows) covers the arrows start .. start + arrows - 1, so
    two consecutive relations share max(0, start + arrows - nextStart) of them.
    """
    relations = lnaMoves.relationsOf(relLengths)
    return tuple(max(0, start + arrows - nextStart)
                 for (start, arrows), (nextStart, _) in zip(relations, relations[1:]))


def maxOverlap(relLengths):
    """The largest overlap between consecutive relations, or 0 if there is none.

    `maxOverlap <= 1` is the almost separate condition, which is why this is the
    coordinate everything below is measured in.
    """
    profile = overlapProfile(relLengths)
    return max(profile) if profile else 0


def overlapRuns(relLengths, threshold = 2):
    """Maximal runs of relations linked by an overlap of at least `threshold`.

    Returned as a list of lists of (start, arrows), longest-to-shortest order
    not imposed -- they come out left to right.  Relations sharing fewer arrows
    than the threshold start a new run, so a run of length one is a relation
    that overlaps nothing heavily.

    The runs are what F-021 is stated in terms of: a run of two is frozen in the
    interior, a run of three or more is not.
    """
    relations = lnaMoves.relationsOf(relLengths)
    if not relations:
        return []
    profile = overlapProfile(relLengths)
    runs = [[relations[0]]]
    for relation, overlap in zip(relations[1:], profile):
        if overlap >= threshold:
            runs[-1].append(relation)
        else:
            runs.append([relation])
    return runs


def isAlmostSeparate(length, relLengths):
    """Whether the quipu theorem names this LNA's class outright."""
    return quipuForms.quipuForAlmostSeparateLNA(length, list(relLengths)) is not None


def moveOrbits(length, rules = None):
    """Partition the LNAs of a length into orbits under the move rules.

    The rules are applied as rewrites on the relation lengths -- no mutation is
    computed, which is legitimate exactly because every rule in the table has
    been verified against the engine.  That makes the whole partition a matter
    of seconds where running the mutations would take hours.

    Returns (list of relation-length tuples, dict from a representative to the
    members of its orbit).
    """
    rules = lnaMoves.ALL_MOVES if rules is None else rules
    lnas = [tuple(algebra.relLengths)
            for algebra in nakayama.LinearNakayamaAlgebra.allOfLength(length)]
    index = {lna: position for position, lna in enumerate(lnas)}
    parent = list(range(len(lnas)))

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for lna in lnas:
        for reached in lnaMoves.rewritesOf(length, list(lna), rules):
            if reached not in index:
                continue
            first, second = find(index[lna]), find(index[reached])
            if first != second:
                parent[first] = second

    orbits = collections.defaultdict(list)
    for lna in lnas:
        orbits[find(index[lna])].append(lna)
    return lnas, orbits


def coverage(length, rules = None):
    """What the quipu theorem plus the move orbits place, broken down by overlap.

    An LNA is *seeded* when the theorem names it, and *covered* when its move
    orbit contains a seeded one -- which is exactly what
    `classification.seedTableFromQuipuTheorem` fills in before any search runs.
    Everything else needs a search.

    Returns a dict with the counts and, under 'byMaxOverlap', a mapping from
    max overlap to (total, covered, uncovered) at that overlap.
    """
    lnas, orbits = moveOrbits(length, rules)
    seeded = {lna for lna in lnas if isAlmostSeparate(length, lna)}
    covered = set()
    for members in orbits.values():
        if any(member in seeded for member in members):
            covered.update(members)
    byMaxOverlap = {}
    for lna in lnas:
        overlap = maxOverlap(lna)
        total, isCovered, _ = byMaxOverlap.get(overlap, (0, 0, 0))
        total += 1
        isCovered += lna in covered
        byMaxOverlap[overlap] = (total, isCovered, total - isCovered)
    return {
        'length': length,
        'lnas': lnas,
        'orbits': orbits,
        'seeded': seeded,
        'covered': covered,
        'uncovered': [lna for lna in lnas if lna not in covered],
        'byMaxOverlap': dict(sorted(byMaxOverlap.items())),
    }


def blockingCores(length, uncovered):
    """The heavily overlapping sub-patterns of the LNAs a search still has to place.

    For each run of relations overlapping in two or more arrows, the run itself,
    normalised so its leftmost relation starts at arrow 1.  Counting these says
    what discovery should be aimed at, rather than at small patterns chosen for
    being cheap to enumerate (research H-003, NOTES idea 19).
    """
    counts = collections.Counter()
    for relLengths in uncovered:
        for run in overlapRuns(list(relLengths)):
            if len(run) < 2:
                continue
            shift = run[0][0] - 1
            counts[tuple((start - shift, arrows) for start, arrows in run)] += 1
    return counts
