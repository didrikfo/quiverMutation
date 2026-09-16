"""Relations of two arrows are free, and what that does to the search.

`corollary:lengthtworelations` of arXiv:2310.08346 says a relation of two arrows
does not affect an algebra's derived equivalence type, so one may be added to or
removed from an LNA anywhere it fits.  That is not a mutation and has no
mutation sequence behind it -- which is exactly why it belongs here and not in
`lnaMoves`, whose every rule is a rewrite the engine has been made to perform.
Everything in this module is a statement about *derived* equivalence; nothing in
it may be used to claim two algebras are in the same mutation class.

Three things follow, and research F-028 is the measurement of all three.

**A relation of two arrows never overlaps a neighbour in more than one arrow.**
Admissibility makes relation starts and ends both strictly increase, so a
relation `(s, 2)` and a later `(t, l)` have `s < t`, giving an overlap of
`s + 2 - t <= 1`; and an earlier `(r, k)` has `r + k < s + 2`, giving `r + k - s
<= 1` the same way.  So `isAlmostSeparate` cannot be broken by a relation of two
arrows, and cannot be repaired by deleting one either: `stripLengthTwo` never
changes whether the quipu theorem names an LNA outright.

**The value of the move is bridging, not naming.**  Since stripping names
nothing on its own, everything it buys is in joining one move orbit to another:
an LNA that no rule reaches strips to one that the rules do reach.

**The reduced space is one vertex smaller.**  Shortening every relation by one
arrow is a bijection from the stripped LNAs of length `n` to *all* the LNAs of
length `n - 1` (`shortenByOneArrow`), so quotienting by the free move divides
the space by the Catalan ratio, about 3.5 by `n = 12` and tending to 4.
"""

import collections

from . import edgeMoves
from . import lnaMoves
from . import nakayama
from . import overlap


def lengthTwoRelations(relLengths):
    """The start vertices of the relations of two arrows."""
    return tuple(start for start, arrows in enumerate(relLengths, start = 1)
                 if arrows == 2)


def stripLengthTwo(relLengths):
    """The same LNA with every relation of two arrows deleted.

    Deleting relations can never break admissibility -- a smaller generating set
    of a monomial ideal stays minimal -- so the result is always an LNA of the
    same length, and it is the canonical representative of its free-move class.
    """
    return tuple(0 if arrows == 2 else arrows for arrows in relLengths)


def isReduced(relLengths):
    """Whether the LNA is its own canonical representative."""
    return 2 not in tuple(relLengths)


def shortenByOneArrow(relLengths):
    """A reduced LNA of length n as an arbitrary LNA of length n - 1.

    Every relation of a reduced LNA has at least three arrows, so taking one
    arrow off each leaves at least two, and the starts and ends still strictly
    increase.  The last entry of a reduced LNA is always 0 -- a relation
    starting at vertex n - 1 would need two arrows and have only one -- so
    dropping it loses nothing, and the map is onto.
    """
    relLengths = tuple(relLengths)
    if not isReduced(relLengths):
        raise ValueError("only a reduced LNA shortens: {0!r}".format(relLengths))
    return tuple(arrows - 1 if arrows else 0 for arrows in relLengths[:-1])


def reducedForms(length):
    """Every canonical representative at a length, in the table's order."""
    seen = set()
    forms = []
    for relLengths in nakayama.allRelationLengths(length):
        reduced = stripLengthTwo(relLengths)
        if reduced not in seen:
            seen.add(reduced)
            forms.append(reduced)
    return forms


def derivedOrbits(length, rules = None, free = True, edges = False):
    """Partition the LNAs of a length under the move rules and what else is asked.

    The same union-find as `overlap.moveOrbits`, with each LNA additionally
    joined to its stripped form when `free`, and to whatever `edgeMoves` reaches
    when `edges`.  With `free` the result is a partition into sets known to be
    derived equivalent, which is coarser than the mutation-class partition that
    `overlap.moveOrbits` gives -- see the module docstring.  With `edges` alone
    it is still a mutation-class partition, since those moves carry sequences.
    """
    rules = lnaMoves.ALL_MOVES if rules is None else rules
    lnas = list(nakayama.allRelationLengths(length))
    index = {lna: position for position, lna in enumerate(lnas)}
    parent = list(range(len(lnas)))

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(first, second):
        first, second = find(first), find(second)
        if first != second:
            parent[first] = second

    for lna in lnas:
        for reached in lnaMoves.rewritesOf(length, list(lna), rules):
            if reached in index:
                union(index[lna], index[reached])
        if free:
            reduced = stripLengthTwo(lna)
            if reduced != lna and reduced in index:
                union(index[lna], index[reduced])
        if edges:
            for reached, _sequence in edgeMoves.rewritesOf(length, lna):
                if reached in index:
                    union(index[lna], index[reached])

    orbits = collections.defaultdict(list)
    for lna in lnas:
        orbits[find(index[lna])].append(lna)
    return lnas, orbits


def coverage(length, rules = None, free = True, edges = False):
    """`overlap.coverage`, with the free move and the edge moves as asked for.

    The same shape of answer, so any two can be printed side by side: what the
    quipu theorem seeds, what the orbits then cover, and what is left for a
    search.
    """
    lnas, orbits = derivedOrbits(length, rules, free, edges)
    seeded = {lna for lna in lnas if overlap.isAlmostSeparate(length, lna)}
    covered = set()
    for members in orbits.values():
        if any(member in seeded for member in members):
            covered.update(members)
    byMaxOverlap = {}
    for lna in lnas:
        value = overlap.maxOverlap(lna)
        total, isCovered, _ = byMaxOverlap.get(value, (0, 0, 0))
        total += 1
        isCovered += lna in covered
        byMaxOverlap[value] = (total, isCovered, total - isCovered)
    return {
        'length': length,
        'lnas': lnas,
        'orbits': orbits,
        'seeded': seeded,
        'covered': covered,
        'uncovered': [lna for lna in lnas if lna not in covered],
        'byMaxOverlap': dict(sorted(byMaxOverlap.items())),
    }
