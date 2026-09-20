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

from . import doubleMutation
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


def derivedOrbits(length, rules = None, free = True, edges = False, doubles = False):
    """Partition the LNAs of a length under the move rules and what else is asked.

    The same union-find as `overlap.moveOrbits`, with each LNA additionally
    joined to its stripped form when `free`, to whatever `edgeMoves` reaches
    when `edges`, and to whatever `doubleMutation` reaches when `doubles`.  With `free` the result is a partition into sets known to be
    derived equivalent, which is coarser than the mutation-class partition that
    `overlap.moveOrbits` gives -- see the module docstring.  With `edges` alone
    it is still a mutation-class partition, since those moves carry sequences,
    and the same holds for `doubles`.  An empty `rules` skips the table, which is
    most of the cost from n = 10 up.
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
        if doubles:
            for reached, _sequence in doubleMutation.rewritesOf(length, lna):
                union(index[lna], index[reached])

    orbits = collections.defaultdict(list)
    for lna in lnas:
        orbits[find(index[lna])].append(lna)
    return lnas, orbits


def coverage(length, rules = None, free = True, edges = False, doubles = False):
    """`overlap.coverage`, with the free move and the edge moves as asked for.

    The same shape of answer, so any two can be printed side by side: what the
    quipu theorem seeds, what the orbits then cover, and what is left for a
    search.
    """
    lnas, orbits = derivedOrbits(length, rules, free, edges, doubles)
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


def movesFrom(length, relLengths, rules = None, free = False, edges = True,
              doubles = True):
    """Every LNA one move away from this one, as a list of relation-length tuples."""
    rules = lnaMoves.ALL_MOVES if rules is None else rules
    current = tuple(relLengths)
    reached = [tuple(name) for name in lnaMoves.rewritesOf(length, list(current), rules)]
    if edges:
        reached += [tuple(name) for name, _sequence in edgeMoves.rewritesOf(length, current)]
    if doubles:
        reached += [tuple(name) for name, _sequence in doubleMutation.rewritesOf(length, current)]
    if free:
        stripped = stripLengthTwo(current)
        if stripped != current:
            reached.append(stripped)
    return reached


def orbitOf(length, relLengths, rules = None, free = False, edges = True,
            doubles = True, limit = 100000, target = None, stopWhen = None):
    """The move orbit of **one** LNA, without partitioning the whole length.

    `derivedOrbits` answers the same question for every LNA at once, which is the
    right shape below `n = 12` and unaffordable above it: at `n = 13` there are
    208012 rows and the interesting configurations are a handful of them.  This
    walks outwards from one row instead, so a single long quiver can be asked
    about directly.

    Returns the set of relation-length tuples reached.  With `free` the orbit is
    of derived equivalence and with it off of mutation equivalence, exactly as in
    `derivedOrbits`.

    **It walks forwards only**, which `derivedOrbits` does not: that one unions
    the rewrites of *every* LNA, so it sees a move into this row from a row this
    row cannot reach.  The moves are not symmetric as rewrites -- `000000` has no
    relation for any of them to act on, while `000002` double-mutates onto it --
    so a reachable target is a proof and an unreachable one is not a disproof.
    Ask it in the direction the moves run: from the LNA that has the relation
    towards the one that does not.

    `target` stops the walk the moment that LNA is reached, which is what a
    membership question wants: the orbit of a long quiver runs to thousands of
    rows and every step costs the whole rule table, where the answer is usually
    in the first few hundred.

    `stopWhen` is the same economy for a membership question with no single
    target: a predicate on a row, and the walk returns as soon as one satisfies
    it.  Most membership questions here are of that kind -- "does this orbit
    hold *any* almost separate LNA" -- and enumerating the orbit to the cap and
    only then looking costs the whole walk to answer a question the first few
    hundred rows had already answered.

    **What the returned set means depends on why the walk stopped.**  Stopped by
    `target` or `stopWhen` it is a prefix of the orbit and its size is how far
    the walk got, not how big the orbit is.  Run to an empty frontier it is the
    whole forward orbit, and *that* is the only case in which a row's absence
    from it says anything.  Callers that need to tell the three apart should use
    `orbitReport`, which says which happened.
    """
    return orbitReport(length, relLengths, rules = rules, free = free,
                       edges = edges, doubles = doubles, limit = limit,
                       target = target, stopWhen = stopWhen).rows


class OrbitWalk(object):
    """The rows a walk reached, and why it stopped -- which is the honest part.

    `rows` is what `orbitOf` returns.  `stoppedBy` is one of:

    * `'closed'` -- the frontier emptied.  `rows` is the entire forward orbit,
      and a row's *absence* from it is a statement about the moves.
    * `'found'`  -- a `target` or a `stopWhen` matched, and `found` is the row
      that did.  `rows` is a prefix of the orbit; absence means nothing.
    * `'cap'`    -- `limit` rows were reached first.  The budget was measured
      and nothing else was; this is the mistake E-037 recorded as a result.

    `orbitOf` returns only `rows`, so a caller that reads a short return as "the
    orbit closed" is right only while nothing stops the walk early.  Anything
    that turns an empty search into a negative answer should ask here instead.
    """

    def __init__(self, rows, stoppedBy, found = None):
        self.rows = rows
        self.stoppedBy = stoppedBy
        self.found = found

    @property
    def closed(self):
        return self.stoppedBy == 'closed'

    def __len__(self):
        return len(self.rows)


def orbitReport(length, relLengths, rules = None, free = False, edges = True,
                doubles = True, limit = 100000, target = None, stopWhen = None):
    """`orbitOf`'s walk, with why it stopped.  See `OrbitWalk`."""
    rules = lnaMoves.ALL_MOVES if rules is None else rules
    start = tuple(relLengths)
    seen = {start}
    wanted = None if target is None else tuple(target)
    if wanted is not None and wanted == start:
        return OrbitWalk(seen, 'found', start)
    if stopWhen is not None and stopWhen(start):
        return OrbitWalk(seen, 'found', start)
    frontier = [start]
    while frontier and len(seen) < limit:
        current = frontier.pop()
        for name in movesFrom(length, current, rules, free, edges, doubles):
            if name not in seen:
                seen.add(name)
                frontier.append(name)
                if wanted is not None and name == wanted:
                    return OrbitWalk(seen, 'found', name)
                if stopWhen is not None and stopWhen(name):
                    return OrbitWalk(seen, 'found', name)
    return OrbitWalk(seen, 'closed' if not frontier else 'cap')


def movesJoin(length, first, second, rules = None, free = False, edges = True,
              doubles = True, limit = 100000):
    """Do the moves join these two LNAs?  Walked from **both** ends at once.

    Every move is a mutation equivalence, and equivalence is symmetric even where
    the rewrite that certifies it is not, so two rows whose forward orbits share
    a row are equivalent.  That is worth having because `orbitOf` walking one way
    answers a membership question only while the orbit stays small enough to
    enumerate, and the orbits of the long quivers are not: at `n = 15` a walk that
    misses its target after twenty thousand rows has proved nothing at all.
    Meeting in the middle reaches the same rows for the square root of the work,
    exactly as `search.meetingPoints` does for quivers.

    Returns the row the two walks meet at, or `None` if neither orbit grew past
    `limit` rows without meeting.  `None` is still not a disproof: it says the
    moves did not join them within the budget.  With `free` on, the deletion of a
    two-arrow relation joins the walk and what is proved is a *derived*
    equivalence rather than a mutation one, exactly as in `derivedOrbits`.
    """
    ends = [tuple(first), tuple(second)]
    seen = [{ends[0]}, {ends[1]}]
    frontier = [[ends[0]], [ends[1]]]
    if ends[0] == ends[1]:
        return ends[0]
    while frontier[0] or frontier[1]:
        side = 0 if len(seen[0]) <= len(seen[1]) else 1
        if not frontier[side]:
            side = 1 - side
        if len(seen[side]) >= limit:
            return None
        current = frontier[side].pop()
        for name in movesFrom(length, current, rules, free, edges, doubles):
            if name in seen[side]:
                continue
            if name in seen[1 - side]:
                return name
            seen[side].add(name)
            frontier[side].append(name)
    return None
