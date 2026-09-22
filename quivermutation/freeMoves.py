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


def addableLengthTwo(length, relLengths):
    """The start vertices at which a relation of two arrows can be added.

    The inverse of the single deletions `stripLengthTwo` composes: an empty
    vertex that a two-arrow relation can start at without breaking the
    strictly increasing starts and ends.
    """
    relLengths = tuple(relLengths)
    found = []
    for position, arrows in enumerate(relLengths):
        if arrows:
            continue
        candidate = list(relLengths)
        candidate[position] = 2
        if lnaMoves.isAdmissible(length, candidate):
            found.append(position + 1)
    return tuple(found)


def reducedMovesFrom(length, relLengths, rules = None, edges = True, doubles = True):
    """Every reduced LNA one step away from this one's reduced form.

    The walk on the quotient by the free move, which is what the free move
    *means*: a relation of two arrows can be deleted **or added**, so an LNA and
    its stripped form are one state, and a state is named by its reduced
    representative.  A step is one mutation move out of the reduced row, or out
    of the reduced row with a single two-arrow relation added, stripped again.

    `movesFrom(free = True)` only ever deletes.  That makes the free move one-way
    in a walk, and the walk then answers differently for two rows that are one
    derived class: at `n = 11` the row `0404…` has **no** move at all while
    `2404…` reaches an almost separate LNA in nine rows (E-049).  This is what
    adding a relation buys -- the two-arrow relation is the spectator or the
    companion a rule needs to fire (F-023, F-032).

    Adding one relation at a time is not the whole class: a row with two
    two-arrow relations is two steps away and is not tried as a starting point
    for a move.  So this is sound -- every step is a derived equivalence -- and
    still a statement about a move set, never a proof of absence.
    """
    rules = lnaMoves.ALL_MOVES if rules is None else rules
    current = stripLengthTwo(relLengths)
    starts = [current]
    for vertex in addableLengthTwo(length, current):
        decorated = list(current)
        decorated[vertex - 1] = 2
        starts.append(tuple(decorated))
    reached = []
    seen = {current}
    for start in starts:
        for name in movesFrom(length, start, rules, free = False, edges = edges,
                              doubles = doubles):
            reduced = stripLengthTwo(name)
            if reduced not in seen:
                seen.add(reduced)
                reached.append(reduced)
    return reached


# `free = REDUCED` in `movesFrom`, `orbitReport` and `movesJoin` walks the
# quotient by the free move: every state is a reduced row, named by
# `stripLengthTwo`.  `free = True` is the one-way deletion it always was.
REDUCED = 'reduced'


def _startOf(relLengths, free):
    return stripLengthTwo(relLengths) if free == REDUCED else tuple(relLengths)


def movesFrom(length, relLengths, rules = None, free = False, edges = True,
              doubles = True):
    """Every LNA one move away from this one, as a list of relation-length tuples."""
    if free == REDUCED:
        return reducedMovesFrom(length, relLengths, rules, edges, doubles)
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
    start = _startOf(relLengths, free)
    seen = {start}
    wanted = None if target is None else _startOf(target, free)
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
    ends = [_startOf(first, free), _startOf(second, free)]
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


# ---------------------------------------------------------------------------
# A walk that remembers what every earlier walk learned
# ---------------------------------------------------------------------------


def mirrorRow(length, relLengths):
    """The relation dual of a row: vertex `v` goes to `n + 1 - v`.

    It keeps the derived class (F-026) and the move set is closed under it, so
    a row and its mirror are one question.
    """
    relLengths = tuple(relLengths)
    mirrored = [0] * len(relLengths)
    for position, arrows in enumerate(relLengths):
        if arrows:
            mirrored[length - (position + 1) - arrows] = arrows
    return tuple(mirrored)


SHARED = 'shared'


class SharedWalk(object):
    """Walks out of many rows of one length, each using what the others found.

    Every move, the free move in both directions and the relation dual is an
    equivalence, so what a walk learns is about a *class*, not a row:

    * **inside spreads both ways.**  Every row a walk passes through is in the
      class of its start, so once one of them is shown to be in a quipu class,
      all of them are, including rows that only reach it and rows only reached
      from it.  A union-find over classes, keyed by the stripped row and its
      mirror, carries that; a walk stops the moment it touches a class already
      known to be inside.
    * **outside is directional and is only cached as such.**  A closed forward
      orbit without an almost separate row is remembered row by row, and a
      later walk that reaches one of its rows does not expand it: everything
      beyond it is known and holds nothing (F-051's cache).

    `verdict` walks in two phases: the plain walk, which is cheap per row, then,
    only if that closes without a certificate, the reduced walk
    (`reducedMovesFrom`), which can add a two-arrow relation partway along.
    At `n = 11` and 12 this gives the reduced walk's verdict on every placement
    of the core census, at a thirtieth of the reduced walk's cost and a tenth of
    the plain one's (E-050).

    An `outside` given early can become `inside` when a later walk joins its
    class to an inside one.  `verdict` returns the earlier starts it promoted,
    so a caller writing verdicts as it goes can record the correction.  The
    answers are therefore at least as strong as a walk from each row alone, and
    how much stronger depends on which rows were walked first.
    """

    def __init__(self, length, limit = 20000):
        self.length = length
        self.limit = limit
        self._parent = {}
        self._inside = set()
        self._closedPlain = set()
        self._closedReduced = set()
        self._outsideStarts = {}

    def classKey(self, relLengths):
        reduced = stripLengthTwo(relLengths)
        return min(reduced, mirrorRow(self.length, reduced))

    def _find(self, key):
        parent = self._parent
        parent.setdefault(key, key)
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def _union(self, first, second):
        first, second = self._find(first), self._find(second)
        if first != second:
            self._parent[first] = second
            if first in self._inside:
                self._inside.discard(first)
                self._inside.add(second)

    def isInside(self, relLengths):
        """Whether this row's class is known to hold an almost separate row."""
        return self._find(self.classKey(relLengths)) in self._inside

    def _markInside(self, relLengths):
        self._inside.add(self._find(self.classKey(relLengths)))

    def _walk(self, start, step, closed):
        seen = {start}
        frontier = [start]
        while frontier and len(seen) < self.limit:
            current = frontier.pop()
            if current in closed:
                continue
            for name in step(current):
                if name in seen:
                    continue
                seen.add(name)
                self._union(self.classKey(current), self.classKey(name))
                if overlap.isAlmostSeparate(self.length, name) or self.isInside(name):
                    return 'found', name, seen
                frontier.append(name)
        return ('closed' if not frontier else 'cap'), None, seen

    def verdict(self, relLengths, label = None):
        """(verdict, how, rows walked, the row that settled it, promoted labels).

        `how` is `'theorem'`, `'shared'` (the class was already known inside),
        `'plain'` or `'reduced'` (the phase that found a certificate),
        `'closed'` or `'cap'`.  `label` names this start in the promotions a
        later call returns; it defaults to the row.
        """
        length = self.length
        start = tuple(relLengths)
        label = start if label is None else label
        if overlap.isAlmostSeparate(length, start):
            self._markInside(start)
            return 'inside', 'theorem', 1, start, self._promotions()
        if self.isInside(start):
            return 'inside', 'shared', 0, None, []
        walked = 0
        phases = (
            (start, lambda row: movesFrom(length, row, free = True), self._closedPlain, 'plain'),
            (stripLengthTwo(start), lambda row: reducedMovesFrom(length, row),
             self._closedReduced, 'reduced'),
        )
        capped = False
        for phaseStart, step, closed, name in phases:
            stopped, found, seen = self._walk(phaseStart, step, closed)
            walked += len(seen)
            if stopped == 'found':
                self._union(self.classKey(start), self.classKey(found))
                self._markInside(start)
                return 'inside', name, walked, found, self._promotions()
            if stopped == 'closed':
                closed |= seen
            else:
                capped = True
        self._outsideStarts[label] = start
        return ('undecided' if capped else 'outside'), \
            ('cap' if capped else 'closed'), walked, None, []

    def markInside(self, relLengths):
        """Record a certificate found some other way (a join), and promote."""
        self._markInside(relLengths)
        return self._promotions()

    def _promotions(self):
        promoted = [label for label, row in self._outsideStarts.items()
                    if self.isInside(row)]
        for label in promoted:
            del self._outsideStarts[label]
        return promoted
