"""Short mutation sequences that take one LNA to another.

The depth-first search is a blunt instrument: it explores every legal mutation to
a fixed depth and reports whatever LNAs it happens to land on.  But many of the
LNA-to-LNA links follow patterns that hold for every length and every position --
a pair of relations sliding along the quiver, say -- and a pattern like that
applies in constant time where the search costs a whole subtree.

This module finds those links and looks for the patterns in them.  `movesFrom`
enumerates the LNAs reachable from one LNA by a short mutation sequence;
`discoverTranslationInvariantMoves` looks across positions and lengths for the
ones that are really a single rule.
"""

import io
import contextlib

import networkx as nx

from . import invariants
from . import lines
from . import mutation
from . import nakayama
from . import pathAlgebra



def _quiet(function, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return function(*args, **kwargs)


def asRelLengths(pathAlg, length):
    """The per-vertex relation lengths of a path algebra that is again a line.

    Returns None if the quiver is not the line on `length` vertices, or carries a
    relation that is not a single zero relation -- an LNA has only those.
    """
    quiver = pathAlg.quiver
    if len(quiver.nodes) != length or len(quiver.edges) != length - 1:
        return None
    if list(nx.simple_cycles(quiver)):
        return None
    if nx.dag_longest_path_length(quiver) != length - 1:
        return None
    relabelled, _ = _quiet(lines.relabelLineAlgebra, pathAlg, {})
    relLengths = [0] * (length - 2)
    for rel in relabelled.rels:
        if len(rel) != 1:
            return None
        start, arrows = rel[0][0], len(rel[0]) - 1
        if start > length - 2 or relLengths[start - 1]:
            return None
        relLengths[start - 1] = arrows
    return relLengths


def movesFrom(lna, maxSteps = 2, allowLeft = True):
    """Every LNA reachable from `lna` by at most `maxSteps` mutations.

    Returns a dict from the reached LNA's class name to the shortest mutation
    sequence that gets there, as the list of signed vertices
    quiverMutationAtVertices takes: positive for a right mutation, negative for a
    left one.  The starting LNA itself is excluded.

    Left mutations are included by default.  The class search only walks right
    mutations, which makes its reachability one-way; a move meant to be used as a
    rewrite rule wants both directions.
    """
    length = lna.length
    start = lines.className(lna.relLengths)
    best = {}

    def walk(pathAlg, steps, history):
        if steps == 0:
            return
        dual = _quiet(pathAlgebra.dualPathAlgebra, pathAlg)
        for vertex in pathAlg.vertices():
            directions = []
            # A mutation is only a tilting mutation where the procedure's
            # admissibility condition holds; applying it anywhere else still
            # computes a quiver, but not a derived equivalent one.  Left mutation
            # at v is right mutation at v of the dual, so that is where its
            # condition is tested.
            if _quiet(mutation.mutationIsPossibleAtVertex, pathAlg, vertex):
                directions.append(vertex)
            if allowLeft and _quiet(mutation.mutationIsPossibleAtVertex, dual, vertex):
                directions.append(-vertex)
            for signed in directions:
                nextAlg = _quiet(mutation.quiverMutationAtVertices,
                                 _copy(pathAlg), [signed])
                if nextAlg is None:
                    continue
                relLengths = asRelLengths(nextAlg, length)
                sequence = history + [signed]
                if relLengths is not None:
                    name = lines.className(relLengths)
                    if name != start and (name not in best or len(sequence) < len(best[name])):
                        best[name] = sequence
                walk(nextAlg, steps - 1, sequence)

    walk(_copy(lna), maxSteps, [])
    return best


def _copy(pathAlg):
    import copy
    duplicate = pathAlgebra.PathAlgebra()
    duplicate.quiver = copy.deepcopy(pathAlg.quiver)
    duplicate.rels = copy.deepcopy(pathAlg.rels)
    return duplicate


def shiftRelLengths(relLengths, offset):
    """Slide every relation `offset` vertices to the right, or None if it falls off."""
    shifted = [0] * len(relLengths)
    for index, arrows in enumerate(relLengths):
        if not arrows:
            continue
        moved = index + offset
        if moved < 0 or moved >= len(relLengths):
            return None
        if moved + arrows > len(relLengths) + 1:
            return None
        if shifted[moved]:
            return None
        shifted[moved] = arrows
    return shifted


def shiftSequence(sequence, offset):
    """Shift the vertices of a mutation sequence, keeping the directions."""
    return [v + offset if v > 0 else v - offset for v in sequence]


# ---------------------------------------------------------------------------
# Verified moves
#
# Each of these is a rewrite on the per-vertex relation lengths, together with
# the mutation sequence that realises it, so the class of an LNA can be expanded
# without running a search.  Every one is checked exhaustively against the
# mutation engine in tests/test_lna_moves.py, over every LNA of every length in
# a range, so the stated conditions are exactly the conditions under which the
# rewrite holds.
# ---------------------------------------------------------------------------


def relationsOf(relLengths):
    """(start vertex, number of arrows) for each relation, left to right."""
    return [(index + 1, arrows) for index, arrows in enumerate(relLengths) if arrows]


def arrowSpan(start, arrows):
    """The arrows a relation covers.  Arrow i runs from vertex i to vertex i+1."""
    return set(range(start, start + arrows))


def maximallyOverlappingPairs(relLengths):
    """Every pair of equal-length relations starting at consecutive vertices.

    These overlap in all but one arrow, which is the most two relations of the
    same length can overlap while still starting at different vertices.  Each is
    returned as (index of the first in relationsOf, start vertex, length).
    """
    relations = relationsOf(relLengths)
    pairs = []
    for index in range(len(relations) - 1):
        (firstStart, firstLength) = relations[index]
        (secondStart, secondLength) = relations[index + 1]
        if firstLength == secondLength and secondStart == firstStart + 1:
            pairs.append((index, firstStart, firstLength))
    return pairs


def _pairIsIsolated(relLengths, index, start, length):
    """Whether no other relation shares an arrow with the pair's span."""
    span = arrowSpan(start, length) | arrowSpan(start + 1, length)
    relations = relationsOf(relLengths)
    for other in range(len(relations)):
        if other in (index, index + 1):
            continue
        if arrowSpan(*relations[other]) & span:
            return False
    return True


def slidePair(length, relLengths, start, direction):
    """Slide a maximally overlapping pair one vertex left or right.

    `direction` is -1 to slide left or +1 to slide right.  Returns
    (new relation lengths, mutation sequence), or None if the move does not
    apply: no such pair at `start`, another relation shares an arrow with it, or
    there is no room to move.

    The mutation sequence is two mutations at one vertex:

    * left  -- two *right* mutations at `start`, the source of the first relation;
    * right -- two *left* mutations at `start + length + 1`, the target of the
      second relation.

    The two are inverse, as they must be, since left mutation at a vertex undoes
    right mutation at it.
    """
    pair = next((p for p in maximallyOverlappingPairs(relLengths) if p[1] == start), None)
    if pair is None:
        return None
    index, start, pairLength = pair
    if not _pairIsIsolated(relLengths, index, start, pairLength):
        return None

    moved = list(relLengths)
    moved[start - 1] = 0
    moved[start] = 0
    if direction < 0:
        if start - 1 < 1:
            return None
        moved[start - 2] = pairLength
        moved[start - 1] = pairLength
        sequence = [start, start]
    else:
        if start + 2 > len(relLengths) or start + 2 + pairLength > length:
            return None
        moved[start] = pairLength
        moved[start + 1] = pairLength
        sequence = [-(start + pairLength + 1), -(start + pairLength + 1)]
    return moved, sequence


def standardise(pathAlg):
    """A copy of a line quiver renumbered to 1 -> 2 -> ... -> n, and the numbering.

    Returns (standardised path algebra, relation lengths, numbering), where
    numbering[p] is the label that the vertex at standard position p carries in
    the quiver as given.  (None, None, None) if the quiver is not a line.

    Note that lines.relabelLineAlgebra renumbers its argument *in place*, so this
    works on a copy.  Getting that wrong is what made composed move sequences
    come out wrong: the numbering reported no longer described the algebra being
    carried forward.

    The renumbering is the subtle part of composing moves at all.
    quiverMutationAtVertex keeps every vertex label -- the mutated vertex i
    becomes i* but is still called i -- so labels never permute.  What changes is
    the *order* in which those labels appear along the line, so a move that wants
    to mutate at "the source of the first relation" has to be told which label
    that vertex has now.
    """
    length = len(pathAlg.quiver.nodes)
    duplicate = _copy(pathAlg)
    relabelled, numbering = _quiet(lines.relabelLineAlgebra, duplicate, {})
    relLengths = [0] * (length - 2)
    for rel in relabelled.rels:
        if len(rel) != 1:
            return None, None, None
        startVertex, arrowCount = rel[0][0], len(rel[0]) - 1
        if startVertex > length - 2 or relLengths[startVertex - 1]:
            return None, None, None
        relLengths[startVertex - 1] = arrowCount
    return relabelled, relLengths, numbering


def closureUnderMoves(length, relLengths, maxIterations = 10000):
    """Every LNA reachable from this one by any number of verified moves.

    All of them are in the same derived equivalence class, with a mutation path
    that proves it.  Returns a dict from class name to (mutation sequence,
    numbering), in the shape the depth-first search produces: the sequence is in
    the *original* vertex labels, and numbering[p] is the label at standard
    position p of the reached line.

    This is the cheap expansion -- a whole orbit of the class for the cost of
    counting relations, where the search would have to walk a tree to find the
    same members.

    Each step works on the standardised LNA, whose own labels are its standard
    positions, and translates the move's vertices into original labels through
    the numbering accumulated so far.
    """
    startName = lines.className(relLengths)
    identity = {position: position for position in range(1, length + 1)}
    results = {startName: ([], identity)}
    frontier = [(nakayama.LinearNakayamaAlgebra(length, relLengths), relLengths, [], identity)]
    iterations = 0
    while frontier and iterations < maxIterations:
        iterations += 1
        standardAlg, current, path, numbering = frontier.pop()
        for name, sequence in movesByRule(length, current).items():
            if name in results:
                continue
            # standardAlg's own labels are the standard positions, so the move
            # applies to it verbatim; the path needs the original labels.
            mutated = _quiet(mutation.quiverMutationAtVertices, _copy(standardAlg), list(sequence))
            nextAlg, nextLengths, stepNumbering = standardise(mutated)
            if nextLengths is None or lines.className(nextLengths) != name:
                continue
            inOriginalLabels = [
                numbering[v] if v > 0 else -numbering[-v] for v in sequence
            ]
            nextNumbering = {
                position: numbering[stepNumbering[position]]
                for position in range(1, length + 1)
            }
            results[name] = (path + inOriginalLabels, nextNumbering)
            frontier.append((nextAlg, nextLengths, path + inOriginalLabels, nextNumbering))
    return results


def orbitOf(length, relLengths):
    """Just the class names reachable by verified moves, as a sorted list."""
    return sorted(closureUnderMoves(length, relLengths))


# ---------------------------------------------------------------------------
# Discovering moves
#
# Rather than guess a rule and test it, enumerate the short mutation sequences
# that take one LNA to another, describe each as a *local* rewrite -- a window of
# the quiver, what the relations inside it become, and where the mutations happen
# relative to the window -- and report the descriptions that recur across
# lengths and positions.  Those are the candidate rules; each still has to be
# verified exhaustively before being trusted.
# ---------------------------------------------------------------------------


def _relationsInWindow(relations, lo, hi):
    """The relations whose arrows fall inside [lo, hi], as offsets from lo."""
    inside = []
    for start, arrows in relations:
        span = arrowSpan(start, arrows)
        if min(span) >= lo and max(span) <= hi:
            inside.append((start - lo, arrows))
    return tuple(sorted(inside))


def describeLink(length, before, after, sequence):
    """Describe a link between two LNAs as a local rewrite, or None.

    The window is the smallest interval of arrows containing every relation that
    differs between the two, and every vertex the sequence mutates at.  A link is
    only described if all relations outside the window are identical and all
    relations touching the window lie entirely inside it -- otherwise the rewrite
    is not local and cannot be stated as a rule.

    Returns (window width, relations before, relations after, sequence offsets),
    with positions given relative to the window's first arrow.
    """
    beforeRelations = relationsOf(before)
    afterRelations = relationsOf(after)
    changed = set(beforeRelations) ^ set(afterRelations)
    if not changed:
        return None
    touched = set()
    for start, arrows in changed:
        touched |= arrowSpan(start, arrows)
    for vertex in sequence:
        touched.add(abs(vertex))
        touched.add(abs(vertex) - 1)
    lo, hi = min(touched), max(touched)
    lo, hi = max(1, lo), min(length - 1, hi)

    # Every relation meeting the window must be contained in it, on both sides,
    # or the rewrite depends on something it does not describe.
    for relations in (beforeRelations, afterRelations):
        for start, arrows in relations:
            span = arrowSpan(start, arrows)
            if span & set(range(lo, hi + 1)) and not (min(span) >= lo and max(span) <= hi):
                return None
    # And everything outside must be untouched.
    outsideBefore = {r for r in beforeRelations if not (arrowSpan(*r) & set(range(lo, hi + 1)))}
    outsideAfter = {r for r in afterRelations if not (arrowSpan(*r) & set(range(lo, hi + 1)))}
    if outsideBefore != outsideAfter:
        return None

    offsets = tuple(v - lo + 1 if v > 0 else v + lo - 1 for v in sequence)
    return (hi - lo + 1,
            _relationsInWindow(beforeRelations, lo, hi),
            _relationsInWindow(afterRelations, lo, hi),
            offsets)


def discoverMoves(lengths, maxSteps = 2, minOccurrences = 3, allowLeft = True,
                  progress = False):
    """Look for local rewrites that recur across lengths and positions.

    Returns a dict from rewrite description to the list of
    (length, class name before, class name after, window start) it was seen at,
    keeping only those seen at `minOccurrences` distinct places.  A description
    seen at many lengths and offsets is a candidate rule; verifyMove then checks
    whether it holds everywhere it applies.
    """
    seen = {}
    for length in lengths:
        for lna in nakayama.LinearNakayamaAlgebra.allOfLength(length):
            if progress:
                print('  {0} {1}'.format(length, lines.className(lna.relLengths)))
            for name, sequence in movesFrom(lna, maxSteps, allowLeft).items():
                after = [int(c) for c in name]
                description = describeLink(length, lna.relLengths, after, sequence)
                if description is None:
                    continue
                seen.setdefault(description, []).append(
                    (length, lines.className(lna.relLengths), name))
    return {d: places for d, places in seen.items() if len(places) >= minOccurrences}


def formatMove(description):
    """A readable one-line form of a rewrite description."""
    width, before, after, offsets = description
    def relations(rels):
        return " ".join("({0}:{1})".format(start, arrows) for start, arrows in rels) or "-"
    return "window {0} arrows: {1}  ->  {2}   via {3}".format(
        width, relations(before), relations(after), list(offsets))


# ---------------------------------------------------------------------------
# Applying and verifying a discovered rewrite
# ---------------------------------------------------------------------------


def isAdmissible(length, relLengths):
    """Whether these relation lengths describe an LNA of that length.

    Each relation needs at least two arrows and must fit, and the starts and the
    ends must both strictly increase -- the conditions of arXiv:2305.06642.
    """
    relations = relationsOf(relLengths)
    if any(arrows < 2 for _, arrows in relations):
        return False
    if any(start + arrows > length for start, arrows in relations):
        return False
    for earlier, later in zip(relations, relations[1:]):
        if earlier[0] + earlier[1] >= later[0] + later[1]:
            return False
    return True


def matchesAt(length, relLengths, description, windowStart):
    """Whether a rewrite's left-hand side sits at this window position.

    The window must fit in the quiver, the relations inside it must be exactly
    the pattern, and no relation outside may reach into it.
    """
    width, before, _after, _offsets = description
    windowEnd = windowStart + width - 1
    if windowStart < 1 or windowEnd > length - 1:
        return False
    window = set(range(windowStart, windowEnd + 1))
    inside = []
    for start, arrows in relationsOf(relLengths):
        span = arrowSpan(start, arrows)
        if not (span & window):
            continue
        if not (min(span) >= windowStart and max(span) <= windowEnd):
            return False            # a relation straddles the window edge
        inside.append((start - windowStart, arrows))
    return tuple(sorted(inside)) == before


def applyAt(length, relLengths, description, windowStart):
    """Apply a rewrite at a window position: (new relation lengths, sequence).

    None if the result would not be an admissible LNA.
    """
    width, before, after, offsets = description
    windowEnd = windowStart + width - 1
    result = list(relLengths)
    for start, arrows in before:
        result[windowStart + start - 1] = 0
    for start, arrows in after:
        position = windowStart + start - 1
        if position < 0 or position >= len(result) or result[position]:
            return None
        result[position] = arrows
    if not isAdmissible(length, result):
        return None
    sequence = [v + windowStart - 1 if v > 0 else v - windowStart + 1 for v in offsets]
    if any(abs(v) < 1 or abs(v) > length for v in sequence):
        return None
    return result, sequence


def isLegalSequence(pathAlg, sequence):
    """Whether every mutation in the sequence is admissible where it is applied.

    The procedure is only a tilting mutation under its admissibility condition;
    outside it, quiverMutationAtVertex still returns a quiver, but not one whose
    path algebra is derived equivalent.  A rewrite built out of illegal steps can
    therefore land on exactly the predicted LNA and still be wrong, which is what
    the Coxeter polynomial check below catches.
    """
    current = _copy(pathAlg)
    for signed in sequence:
        target = current if signed > 0 else _quiet(pathAlgebra.dualPathAlgebra, current)
        if not _quiet(mutation.mutationIsPossibleAtVertex, target, abs(signed)):
            return False
        current = _quiet(mutation.quiverMutationAtVertices, current, [signed])
    return True


def verifyMove(description, lengths, checkCoxeter = True):
    """Check a rewrite against the mutation engine everywhere it applies.

    Three things have to hold, and all three are needed:

    1. the result is the predicted LNA;
    2. every mutation in the sequence is admissible where it lands, so the
       sequence really is a chain of tilting mutations;
    3. the Coxeter polynomial does not move, since it is a derived invariant.

    Checking only (1) admits rules that are simply false: a sequence containing
    an inadmissible step can still produce the predicted relation lengths, and
    then the two algebras are not derived equivalent at all.  When this check was
    first written without (2) and (3), 38 rules passed and 6561 of the 8388 orbit
    members they generated had the wrong Coxeter polynomial.

    Returns (confirmed, failures), each failure being
    (length, before, predicted, actual, reason).
    """
    confirmed = 0
    failures = []
    for length in lengths:
        # Enumerate the admissible LNAs directly.  Filtering the full product of
        # relation lengths instead means 10^8 tuples at length 10, against the
        # 4862 LNAs that actually exist there.
        for algebra in nakayama.LinearNakayamaAlgebra.allOfLength(length):
            relLengths = algebra.relLengths
            for windowStart in range(1, length):
                if not matchesAt(length, relLengths, description, windowStart):
                    continue
                applied = applyAt(length, relLengths, description, windowStart)
                if applied is None:
                    continue
                predicted, sequence = applied
                startAlg = nakayama.LinearNakayamaAlgebra(length, relLengths)
                if not isLegalSequence(startAlg, sequence):
                    failures.append((length, lines.className(relLengths),
                                     lines.className(predicted), None, 'illegal mutation'))
                    continue
                mutated = _quiet(mutation.quiverMutationAtVertices, _copy(startAlg), list(sequence))
                actual = asRelLengths(mutated, length)
                if actual != predicted:
                    failures.append((length, lines.className(relLengths),
                                     lines.className(predicted),
                                     lines.className(actual) if actual else None, 'wrong result'))
                    continue
                if checkCoxeter and not _sameCoxeter(length, relLengths, predicted):
                    failures.append((length, lines.className(relLengths),
                                     lines.className(predicted), lines.className(actual),
                                     'Coxeter polynomial moved'))
                    continue
                confirmed += 1
    return confirmed, failures


def _sameCoxeter(length, before, after):
    import sympy
    first = _quiet(invariants.coxeterPoly, nakayama.LinearNakayamaAlgebra(length, before)).as_expr()
    second = _quiet(invariants.coxeterPoly, nakayama.LinearNakayamaAlgebra(length, after)).as_expr()
    return sympy.expand(first) == sympy.expand(second)


# ---------------------------------------------------------------------------
# The verified rule table
#
# Each entry is (window width in arrows, relations before, relations after,
# mutation sequence), with relation starts and mutation vertices given relative
# to the window's first arrow, and a negative vertex meaning a left mutation.
#
# These were found by discoverMoves over lengths 6 and 7 and then checked by
# verifyMove against the mutation engine over every LNA of lengths 5 to 8 at
# every window position where they apply -- zero failures each.  The count in
# each comment is how many applications were confirmed.
#
# Two of them are the moves that prompted this work:
#
#   (5, ((0,3),(1,3)), ((1,3),(2,3)), (-5,-5))
#       A maximally overlapping pair of equal-length relations slides one arrow
#       right under two left mutations at the second relation's target; the
#       mirror entry slides it left under two right mutations at the first
#       relation's source.
#
#   (5, ((0,3),(1,3),(3,2)), ((0,2),(1,3),(2,3)), (2,2))
#       A relation stays put while the meeting point of the two relations around
#       it moves one arrow left, under two right mutations.
#
# Rules that look right but are not: "a relation of length 2 sitting alone may be
# deleted" holds 63 times and fails 130 times when stated on a two-arrow window.
# It needs a three-arrow window to be true, which is why nothing goes in this
# table without verifyMove.
# ---------------------------------------------------------------------------

VERIFIED_MOVES = [
    (3, ((0, 2),), ((1, 2),), (-3,)),   # 63 confirmed: window 3 arrows: (0:2)  ->  (1:2)   via [-3]
    (3, ((1, 2),), ((0, 2),), (2,)),   # 63 confirmed: window 3 arrows: (1:2)  ->  (0:2)   via [2]
    (4, ((0, 2),), ((2, 2),), (-3, -4)),   # 22 confirmed: window 4 arrows: (0:2)  ->  (2:2)   via [-3, -4]
    (4, ((2, 2),), ((0, 2),), (3, 2)),   # 22 confirmed: window 4 arrows: (2:2)  ->  (0:2)   via [3, 2]
    (5, ((0, 2), (1, 3)), ((0, 3), (1, 3), (2, 3)), (-5, -5)),   # 8 confirmed: window 5 arrows: (0:2) (1:3)  ->  (0:3) (1:3) (2:3)   via [-5, -5]
    (5, ((0, 2), (1, 3), (2, 3)), ((0, 3), (1, 3), (3, 2)), (-5, -3)),   # 8 confirmed: window 5 arrows: (0:2) (1:3) (2:3)  ->  (0:3) (1:3) (3:2)   via [-5, -3]
    (5, ((0, 2), (2, 2)), ((1, 2), (3, 2)), (-3, -5)),   # 8 confirmed: window 5 arrows: (0:2) (2:2)  ->  (1:2) (3:2)   via [-3, -5]
    (5, ((0, 2), (3, 2)), ((1, 2), (2, 2)), (-3, 4)),   # 8 confirmed: window 5 arrows: (0:2) (3:2)  ->  (1:2) (2:2)   via [-3, 4]
    (5, ((0, 3), (1, 3)), ((1, 3), (2, 3)), (-5, -5)),   # 8 confirmed: window 5 arrows: (0:3) (1:3)  ->  (1:3) (2:3)   via [-5, -5]
    (5, ((0, 3), (1, 3), (2, 3)), ((0, 2), (1, 3)), (2, 2)),   # 8 confirmed: window 5 arrows: (0:3) (1:3) (2:3)  ->  (0:2) (1:3)   via [2, 2]
    (5, ((0, 3), (1, 3), (2, 3)), ((1, 3), (3, 2)), (-5, -5)),   # 8 confirmed: window 5 arrows: (0:3) (1:3) (2:3)  ->  (1:3) (3:2)   via [-5, -5]
    (5, ((0, 3), (1, 3), (3, 2)), ((0, 2), (1, 3), (2, 3)), (2, 2)),   # 8 confirmed: window 5 arrows: (0:3) (1:3) (3:2)  ->  (0:2) (1:3) (2:3)   via [2, 2]
    (5, ((1, 2), (2, 2)), ((0, 2), (3, 2)), (2, -5)),   # 8 confirmed: window 5 arrows: (1:2) (2:2)  ->  (0:2) (3:2)   via [2, -5]
    (5, ((1, 2), (3, 2)), ((0, 2), (2, 2)), (2, 4)),   # 8 confirmed: window 5 arrows: (1:2) (3:2)  ->  (0:2) (2:2)   via [2, 4]
    (5, ((1, 3), (2, 3)), ((0, 3), (1, 3)), (2, 2)),   # 8 confirmed: window 5 arrows: (1:3) (2:3)  ->  (0:3) (1:3)   via [2, 2]
    (5, ((1, 3), (3, 2)), ((0, 3), (1, 3), (2, 3)), (2, 2)),   # 8 confirmed: window 5 arrows: (1:3) (3:2)  ->  (0:3) (1:3) (2:3)   via [2, 2]
    (5, ((0, 3), (1, 4)), ((0, 4), (2, 3)), (2, -4, -5)),   # 8 confirmed: window 5 arrows: (0:3) (1:4)  ->  (0:4) (2:3)   via [2, -4, -5]
    (5, ((0, 4), (2, 3)), ((0, 3), (1, 4)), (3, 2, -5)),   # 8 confirmed: window 5 arrows: (0:4) (2:3)  ->  (0:3) (1:4)   via [3, 2, -5]

    # Found by discovery at length 8, where a six-arrow window has room to sit
    # away from both ends, and re-verified over lengths 7 to 10.
    (6, ((0, 2), (1, 3), (2, 4)), ((0, 3), (1, 3), (2, 3), (3, 3)), (-5, -5)),   # 22 confirmed: window 6 arrows: (0:2) (1:3) (2:4)  ->  (0:3) (1:3) (2:3) (3:3)   via [-5, -5]
    (6, ((0, 2), (1, 4)), ((0, 3), (1, 4), (2, 4)), (-6, -6)),   # 22 confirmed: window 6 arrows: (0:2) (1:4)  ->  (0:3) (1:4) (2:4)   via [-6, -6]
    (6, ((0, 2), (1, 4), (2, 4)), ((0, 3), (1, 4), (3, 3)), (-6, -6)),   # 22 confirmed: window 6 arrows: (0:2) (1:4) (2:4)  ->  (0:3) (1:4) (3:3)   via [-6, -6]
    (6, ((0, 2), (1, 4), (3, 3)), ((0, 3), (1, 4), (4, 2)), (-6, -6)),   # 22 confirmed: window 6 arrows: (0:2) (1:4) (3:3)  ->  (0:3) (1:4) (4:2)   via [-6, -6]
    (6, ((0, 2), (2, 2), (3, 2)), ((1, 2), (2, 2), (4, 2)), (-3, -6)),   # 22 confirmed: window 6 arrows: (0:2) (2:2) (3:2)  ->  (1:2) (2:2) (4:2)   via [-3, -6]
    (6, ((0, 2), (2, 2), (4, 2)), ((1, 2), (2, 2), (3, 2)), (-3, 5)),   # 22 confirmed: window 6 arrows: (0:2) (2:2) (4:2)  ->  (1:2) (2:2) (3:2)   via [-3, 5]
    (6, ((0, 2), (3, 2)), ((1, 2), (4, 2)), (-3, -6)),   # 22 confirmed: window 6 arrows: (0:2) (3:2)  ->  (1:2) (4:2)   via [-3, -6]
    (6, ((0, 2), (4, 2)), ((1, 2), (3, 2)), (-3, 5)),   # 22 confirmed: window 6 arrows: (0:2) (4:2)  ->  (1:2) (3:2)   via [-3, 5]
    (6, ((0, 3), (1, 3), (2, 3)), ((0, 4), (2, 3), (3, 3)), (-6, -6)),   # 22 confirmed: window 6 arrows: (0:3) (1:3) (2:3)  ->  (0:4) (2:3) (3:3)   via [-6, -6]
    (6, ((0, 3), (1, 3), (2, 3), (3, 3)), ((0, 2), (1, 3), (2, 4)), (2, 2)),   # 22 confirmed: window 6 arrows: (0:3) (1:3) (2:3) (3:3)  ->  (0:2) (1:3) (2:4)   via [2, 2]
    (6, ((0, 3), (1, 3), (2, 3), (3, 3)), ((0, 4), (2, 3), (4, 2)), (-6, -6)),   # 22 confirmed: window 6 arrows: (0:3) (1:3) (2:3) (3:3)  ->  (0:4) (2:3) (4:2)   via [-6, -6]
    (6, ((0, 3), (1, 3), (2, 4)), ((1, 3), (2, 3), (3, 3)), (-5, -5)),   # 22 confirmed: window 6 arrows: (0:3) (1:3) (2:4)  ->  (1:3) (2:3) (3:3)   via [-5, -5]
    (6, ((0, 3), (1, 4)), ((0, 4), (1, 4), (2, 4)), (-6, -6)),   # 22 confirmed: window 6 arrows: (0:3) (1:4)  ->  (0:4) (1:4) (2:4)   via [-6, -6]
    (6, ((0, 3), (1, 4), (2, 4)), ((0, 2), (1, 4)), (2, 2)),   # 22 confirmed: window 6 arrows: (0:3) (1:4) (2:4)  ->  (0:2) (1:4)   via [2, 2]
    (6, ((0, 3), (1, 4), (2, 4)), ((0, 4), (1, 4), (3, 3)), (-6, -6)),   # 22 confirmed: window 6 arrows: (0:3) (1:4) (2:4)  ->  (0:4) (1:4) (3:3)   via [-6, -6]
    (6, ((0, 3), (1, 4), (3, 3)), ((0, 2), (1, 4), (2, 4)), (2, 2)),   # 22 confirmed: window 6 arrows: (0:3) (1:4) (3:3)  ->  (0:2) (1:4) (2:4)   via [2, 2]
    (6, ((0, 3), (1, 4), (3, 3)), ((0, 4), (1, 4), (4, 2)), (-6, -6)),   # 22 confirmed: window 6 arrows: (0:3) (1:4) (3:3)  ->  (0:4) (1:4) (4:2)   via [-6, -6]
    (6, ((0, 3), (1, 4), (4, 2)), ((0, 2), (1, 4), (3, 3)), (2, 2)),   # 22 confirmed: window 6 arrows: (0:3) (1:4) (4:2)  ->  (0:2) (1:4) (3:3)   via [2, 2]
    (6, ((0, 4), (1, 4)), ((1, 4), (2, 4)), (-6, -6)),   # 22 confirmed: window 6 arrows: (0:4) (1:4)  ->  (1:4) (2:4)   via [-6, -6]
    (6, ((0, 4), (1, 4), (2, 4)), ((0, 3), (1, 4)), (2, 2)),   # 22 confirmed: window 6 arrows: (0:4) (1:4) (2:4)  ->  (0:3) (1:4)   via [2, 2]
    (6, ((0, 4), (1, 4), (2, 4)), ((1, 4), (3, 3)), (-6, -6)),   # 22 confirmed: window 6 arrows: (0:4) (1:4) (2:4)  ->  (1:4) (3:3)   via [-6, -6]
    (6, ((0, 4), (1, 4), (3, 3)), ((0, 3), (1, 4), (2, 4)), (2, 2)),   # 22 confirmed: window 6 arrows: (0:4) (1:4) (3:3)  ->  (0:3) (1:4) (2:4)   via [2, 2]
    (6, ((0, 4), (1, 4), (3, 3)), ((1, 4), (4, 2)), (-6, -6)),   # 22 confirmed: window 6 arrows: (0:4) (1:4) (3:3)  ->  (1:4) (4:2)   via [-6, -6]
    (6, ((0, 4), (1, 4), (4, 2)), ((0, 3), (1, 4), (3, 3)), (2, 2)),   # 22 confirmed: window 6 arrows: (0:4) (1:4) (4:2)  ->  (0:3) (1:4) (3:3)   via [2, 2]
    (6, ((0, 4), (2, 3), (3, 3)), ((0, 3), (1, 3), (2, 3)), (3, 3)),   # 22 confirmed: window 6 arrows: (0:4) (2:3) (3:3)  ->  (0:3) (1:3) (2:3)   via [3, 3]
    (6, ((0, 4), (2, 3), (4, 2)), ((0, 3), (1, 3), (2, 3), (3, 3)), (3, 3)),   # 22 confirmed: window 6 arrows: (0:4) (2:3) (4:2)  ->  (0:3) (1:3) (2:3) (3:3)   via [3, 3]
    (6, ((1, 2), (2, 2), (3, 2)), ((0, 2), (2, 2), (4, 2)), (2, -6)),   # 22 confirmed: window 6 arrows: (1:2) (2:2) (3:2)  ->  (0:2) (2:2) (4:2)   via [2, -6]
    (6, ((1, 2), (2, 2), (4, 2)), ((0, 2), (2, 2), (3, 2)), (2, 5)),   # 22 confirmed: window 6 arrows: (1:2) (2:2) (4:2)  ->  (0:2) (2:2) (3:2)   via [2, 5]
    (6, ((1, 2), (3, 2)), ((0, 2), (4, 2)), (2, -6)),   # 22 confirmed: window 6 arrows: (1:2) (3:2)  ->  (0:2) (4:2)   via [2, -6]
    (6, ((1, 2), (4, 2)), ((0, 2), (3, 2)), (2, 5)),   # 22 confirmed: window 6 arrows: (1:2) (4:2)  ->  (0:2) (3:2)   via [2, 5]
    (6, ((1, 3), (2, 3), (3, 3)), ((0, 3), (1, 3), (2, 4)), (2, 2)),   # 22 confirmed: window 6 arrows: (1:3) (2:3) (3:3)  ->  (0:3) (1:3) (2:4)   via [2, 2]
    (6, ((1, 4), (2, 4)), ((0, 4), (1, 4)), (2, 2)),   # 22 confirmed: window 6 arrows: (1:4) (2:4)  ->  (0:4) (1:4)   via [2, 2]
    (6, ((1, 4), (3, 3)), ((0, 4), (1, 4), (2, 4)), (2, 2)),   # 22 confirmed: window 6 arrows: (1:4) (3:3)  ->  (0:4) (1:4) (2:4)   via [2, 2]
    (6, ((1, 4), (4, 2)), ((0, 4), (1, 4), (3, 3)), (2, 2)),   # 22 confirmed: window 6 arrows: (1:4) (4:2)  ->  (0:4) (1:4) (3:3)   via [2, 2]
    # Interior discovery at three mutations, E-011: patterns planted in the
    # middle of A_13 and A_14 and mutated only nearby, so that no end of the
    # quiver is in reach (H-007).  Each was then verified at
    # `width + 1 .. width + 4` -- the lengths have to follow the window, and
    # verifying these at a fixed 7 to 10 admitted 30 rules of window 9 that are
    # false at length 11 (R-009).  The window-7 and window-8 entries were
    # checked at lengths 11 and 12 as well.
    (5, ((0, 2),), ((3, 2),), (-3, -4, -5)),   # 22 confirmed: window 5 arrows: (0:2)  ->  (3:2)   via [-3, -4, -5]
    (5, ((3, 2),), ((0, 2),), (4, 3, 2)),   # 22 confirmed: window 5 arrows: (3:2)  ->  (0:2)   via [4, 3, 2]
    (6, ((0, 2), (1, 3), (2, 3)), ((0, 3), (1, 3), (4, 2)), (-5, -3, -6)),   # 22 confirmed: window 6 arrows: (0:2) (1:3) (2:3)  ->  (0:3) (1:3) (4:2)   via [-5, -3, -6]
    (6, ((0, 2), (2, 2)), ((1, 2), (4, 2)), (-3, -5, -6)),   # 22 confirmed: window 6 arrows: (0:2) (2:2)  ->  (1:2) (4:2)   via [-3, -5, -6]
    (6, ((0, 2), (2, 3)), ((0, 3), (1, 3), (2, 4)), (-6, 3, -6)),   # 22 confirmed: window 6 arrows: (0:2) (2:3)  ->  (0:3) (1:3) (2:4)   via [-6, 3, -6]
    (6, ((0, 2), (2, 3)), ((1, 3), (2, 3), (3, 3)), (-3, -6, -6)),   # 22 confirmed: window 6 arrows: (0:2) (2:3)  ->  (1:3) (2:3) (3:3)   via [-3, -6, -6]
    (6, ((0, 2), (3, 2)), ((2, 2), (4, 2)), (-3, -4, -6)),   # 22 confirmed: window 6 arrows: (0:2) (3:2)  ->  (2:2) (4:2)   via [-3, -4, -6]
    (6, ((0, 3), (1, 3), (2, 3)), ((1, 3), (4, 2)), (-5, -5, -6)),   # 22 confirmed: window 6 arrows: (0:3) (1:3) (2:3)  ->  (1:3) (4:2)   via [-5, -5, -6]
    (6, ((0, 4), (2, 3)), ((1, 3), (2, 4)), (-6, 3, -6)),   # 22 confirmed: window 6 arrows: (0:4) (2:3)  ->  (1:3) (2:4)   via [-6, 3, -6]
    (6, ((1, 2), (2, 2)), ((0, 2), (4, 2)), (2, -5, -6)),   # 22 confirmed: window 6 arrows: (1:2) (2:2)  ->  (0:2) (4:2)   via [2, -5, -6]
    (6, ((1, 2), (4, 2)), ((0, 2), (2, 2)), (2, 5, 4)),   # 22 confirmed: window 6 arrows: (1:2) (4:2)  ->  (0:2) (2:2)   via [2, 5, 4]
    (6, ((1, 3), (2, 3), (3, 3)), ((0, 2), (2, 3)), (3, 2, 3)),   # 22 confirmed: window 6 arrows: (1:3) (2:3) (3:3)  ->  (0:2) (2:3)   via [3, 2, 3]
    (6, ((1, 3), (2, 3), (4, 2)), ((0, 2), (2, 3), (3, 3)), (3, 2, 3)),   # 22 confirmed: window 6 arrows: (1:3) (2:3) (4:2)  ->  (0:2) (2:3) (3:3)   via [3, 2, 3]
    (6, ((1, 3), (2, 4)), ((0, 4), (2, 3)), (2, -5, 2)),   # 22 confirmed: window 6 arrows: (1:3) (2:4)  ->  (0:4) (2:3)   via [2, -5, 2]
    (6, ((1, 3), (4, 2)), ((0, 3), (1, 3), (2, 3)), (2, 5, 2)),   # 22 confirmed: window 6 arrows: (1:3) (4:2)  ->  (0:3) (1:3) (2:3)   via [2, 5, 2]
    (6, ((1, 3), (4, 2)), ((0, 4), (2, 3), (3, 3)), (2, -5, 2)),   # 22 confirmed: window 6 arrows: (1:3) (4:2)  ->  (0:4) (2:3) (3:3)   via [2, -5, 2]
    (6, ((2, 2), (3, 2)), ((0, 2), (4, 2)), (3, 2, -6)),   # 22 confirmed: window 6 arrows: (2:2) (3:2)  ->  (0:2) (4:2)   via [3, 2, -6]
    (6, ((2, 2), (4, 2)), ((0, 2), (3, 2)), (3, 2, 5)),   # 22 confirmed: window 6 arrows: (2:2) (4:2)  ->  (0:2) (3:2)   via [3, 2, 5]
    (7, ((0, 2), (2, 2), (3, 2)), ((1, 2), (2, 2), (5, 2)), (-3, -6, -7)),   # 3 confirmed: window 7 arrows: (0:2) (2:2) (3:2)  ->  (1:2) (2:2) (5:2)   via [-3, -6, -7]
    (7, ((0, 2), (3, 2)), ((1, 2), (5, 2)), (-3, -6, -7)),   # 3 confirmed: window 7 arrows: (0:2) (3:2)  ->  (1:2) (5:2)   via [-3, -6, -7]
    (7, ((1, 2), (2, 2), (3, 2)), ((0, 2), (2, 2), (5, 2)), (2, -6, -7)),   # 3 confirmed: window 7 arrows: (1:2) (2:2) (3:2)  ->  (0:2) (2:2) (5:2)   via [2, -6, -7]
    (7, ((1, 2), (2, 2), (3, 3)), ((0, 2), (2, 3), (3, 3), (4, 3)), (2, -7, -7)),   # 3 confirmed: window 7 arrows: (1:2) (2:2) (3:3)  ->  (0:2) (2:3) (3:3) (4:3)   via [2, -7, -7]
    (7, ((1, 2), (2, 2), (4, 2)), ((0, 2), (2, 2), (5, 2)), (2, -7)),   # 3 confirmed: window 7 arrows: (1:2) (2:2) (4:2)  ->  (0:2) (2:2) (5:2)   via [2, -7]
    (7, ((1, 2), (2, 2), (4, 2)), ((0, 2), (3, 2), (5, 2)), (2, -5, -7)),   # 3 confirmed: window 7 arrows: (1:2) (2:2) (4:2)  ->  (0:2) (3:2) (5:2)   via [2, -5, -7]
    (7, ((1, 2), (2, 3), (3, 3)), ((0, 2), (3, 3), (4, 3)), (2, -7, -7)),   # 3 confirmed: window 7 arrows: (1:2) (2:3) (3:3)  ->  (0:2) (3:3) (4:3)   via [2, -7, -7]
    (7, ((1, 2), (2, 3), (4, 2)), ((0, 2), (2, 3), (5, 2)), (2, -7)),   # 3 confirmed: window 7 arrows: (1:2) (2:3) (4:2)  ->  (0:2) (2:3) (5:2)   via [2, -7]
    (7, ((1, 2), (3, 2)), ((0, 2), (5, 2)), (2, -6, -7)),   # 3 confirmed: window 7 arrows: (1:2) (3:2)  ->  (0:2) (5:2)   via [2, -6, -7]
    (7, ((1, 2), (3, 2), (4, 2)), ((0, 2), (2, 2), (5, 2)), (2, 4, -7)),   # 3 confirmed: window 7 arrows: (1:2) (3:2) (4:2)  ->  (0:2) (2:2) (5:2)   via [2, 4, -7]
    (7, ((1, 2), (3, 2), (4, 2)), ((0, 2), (3, 2), (5, 2)), (2, -7)),   # 3 confirmed: window 7 arrows: (1:2) (3:2) (4:2)  ->  (0:2) (3:2) (5:2)   via [2, -7]
    (7, ((1, 2), (4, 2)), ((0, 2), (5, 2)), (2, -7)),   # 3 confirmed: window 7 arrows: (1:2) (4:2)  ->  (0:2) (5:2)   via [2, -7]
    (7, ((1, 3), (2, 3), (4, 2)), ((0, 3), (1, 3), (5, 2)), (2, 2, -7)),   # 3 confirmed: window 7 arrows: (1:3) (2:3) (4:2)  ->  (0:3) (1:3) (5:2)   via [2, 2, -7]
    (7, ((1, 3), (3, 2), (4, 2)), ((0, 3), (1, 3), (2, 3), (5, 2)), (2, 2, -7)),   # 3 confirmed: window 7 arrows: (1:3) (3:2) (4:2)  ->  (0:3) (1:3) (2:3) (5:2)   via [2, 2, -7]
    (7, ((2, 2), (3, 2), (4, 2)), ((0, 2), (3, 2), (5, 2)), (3, 2, -7)),   # 3 confirmed: window 7 arrows: (2:2) (3:2) (4:2)  ->  (0:2) (3:2) (5:2)   via [3, 2, -7]
    (7, ((2, 2), (3, 2), (5, 2)), ((0, 2), (3, 2), (4, 2)), (3, 2, 6)),   # 3 confirmed: window 7 arrows: (2:2) (3:2) (5:2)  ->  (0:2) (3:2) (4:2)   via [3, 2, 6]
    (7, ((2, 2), (4, 2)), ((0, 2), (5, 2)), (3, 2, -7)),   # 3 confirmed: window 7 arrows: (2:2) (4:2)  ->  (0:2) (5:2)   via [3, 2, -7]
    (7, ((2, 2), (5, 2)), ((0, 2), (4, 2)), (3, 2, 6)),   # 3 confirmed: window 7 arrows: (2:2) (5:2)  ->  (0:2) (4:2)   via [3, 2, 6]
    (8, ((1, 2), (2, 2), (4, 2)), ((0, 2), (2, 2), (6, 2)), (2, -7, -8)),   # 3 confirmed: window 8 arrows: (1:2) (2:2) (4:2)  ->  (0:2) (2:2) (6:2)   via [2, -7, -8]
    (8, ((1, 2), (2, 3), (4, 2)), ((0, 2), (2, 3), (6, 2)), (2, -7, -8)),   # 3 confirmed: window 8 arrows: (1:2) (2:3) (4:2)  ->  (0:2) (2:3) (6:2)   via [2, -7, -8]
    (8, ((1, 2), (3, 2), (4, 2)), ((0, 2), (3, 2), (6, 2)), (2, -7, -8)),   # 3 confirmed: window 8 arrows: (1:2) (3:2) (4:2)  ->  (0:2) (3:2) (6:2)   via [2, -7, -8]
    (8, ((1, 2), (4, 2)), ((0, 2), (6, 2)), (2, -7, -8)),   # 3 confirmed: window 8 arrows: (1:2) (4:2)  ->  (0:2) (6:2)   via [2, -7, -8]
    (8, ((2, 2), (3, 2), (5, 2)), ((0, 2), (3, 2), (6, 2)), (3, 2, -8)),   # 3 confirmed: window 8 arrows: (2:2) (3:2) (5:2)  ->  (0:2) (3:2) (6:2)   via [3, 2, -8]
    (8, ((2, 2), (3, 3), (5, 2)), ((0, 2), (3, 3), (6, 2)), (3, 2, -8)),   # 3 confirmed: window 8 arrows: (2:2) (3:3) (5:2)  ->  (0:2) (3:3) (6:2)   via [3, 2, -8]
    (8, ((2, 2), (4, 2), (5, 2)), ((0, 2), (4, 2), (6, 2)), (3, 2, -8)),   # 3 confirmed: window 8 arrows: (2:2) (4:2) (5:2)  ->  (0:2) (4:2) (6:2)   via [3, 2, -8]
    (8, ((2, 2), (5, 2)), ((0, 2), (6, 2)), (3, 2, -8)),   # 3 confirmed: window 8 arrows: (2:2) (5:2)  ->  (0:2) (6:2)   via [3, 2, -8]
]


def movesByRule(length, relLengths):
    """Every LNA reachable from this one by a single verified move.

    Returns a dict from the reached class name to the mutation sequence, in the
    same shape as movesFrom, but computed by table lookup rather than by walking
    a mutation tree.
    """
    reached = {}
    for description in VERIFIED_MOVES:
        for windowStart in range(1, length):
            if not matchesAt(length, relLengths, description, windowStart):
                continue
            applied = applyAt(length, relLengths, description, windowStart)
            if applied is None:
                continue
            moved, sequence = applied
            name = lines.className(moved)
            if name != lines.className(relLengths) and name not in reached:
                reached[name] = sequence
    return reached


# ---------------------------------------------------------------------------
# Discovering moves in the interior of a long quiver
#
# The first discovery pass looked at every LNA of lengths 6 and 7, which is
# exactly the wrong place to look: on a quiver that short every vertex is within
# a step or two of an end, so the special cases that apply near the boundary
# apply almost everywhere, and a rule that is really about the interior cannot be
# told apart from one that depends on an end being close by.
#
# The fix is to embed a small pattern of relations in the *middle* of a long
# quiver, with several arrows of empty quiver on each side, and to allow
# mutations only at vertices near the pattern.  That does two things at once: the
# rewrite discovered is genuinely local and position-independent by construction,
# and the branching factor stops depending on the length of the quiver, so
# sequences of four or five mutations become affordable where enumerating over the
# whole quiver would not.
# ---------------------------------------------------------------------------


def embedPattern(length, pattern, offset):
    """Place a pattern of (relative start, arrows) at `offset` arrows in.

    Returns the relation lengths, or None if it does not fit or is inadmissible.
    """
    relLengths = [0] * (length - 2)
    for start, arrows in pattern:
        position = offset + start - 1
        if position < 0 or position >= len(relLengths):
            return None
        if relLengths[position]:
            return None
        relLengths[position] = arrows
    if not isAdmissible(length, relLengths):
        return None
    return relLengths


def patternWidth(pattern):
    """How many arrows a pattern spans."""
    covered = set()
    for start, arrows in pattern:
        covered |= arrowSpan(start, arrows)
    return max(covered) - min(covered) + 1 if covered else 0


def _stateKey(pathAlg):
    arrows = tuple(sorted((a[0], a[1]) for a in pathAlg.quiver.edges))
    rels = tuple(sorted(tuple(tuple(p) for p in sorted(rel)) for rel in pathAlg.rels))
    return (arrows, rels)


def localMutationSequences(length, relLengths, centreLo, centreHi, maxSteps, margin):
    """Mutation sequences near a region, and the LNAs they reach.

    Only vertices within `margin` of [centreLo, centreHi] are mutated, and only
    where the mutation is admissible.  Intermediate quivers already seen are not
    re-explored, which is what keeps four- and five-step sequences affordable.

    Returns a dict from reached relation lengths (as a class name) to the
    shortest sequence found.
    """
    allowed = [v for v in range(max(1, centreLo - margin),
                                min(length, centreHi + 1 + margin) + 1)]
    startName = lines.className(relLengths)
    best = {}
    # Maps an intermediate quiver to the most steps that were still available
    # when it was last explored.  Pruning on mere membership loses paths: a state
    # first reached deep in one branch would block a later branch that reaches it
    # with more steps left, so the search would find *fewer* LNAs at a higher
    # maxSteps than at a lower one.
    seen = {}

    def walk(pathAlg, steps, history):
        if steps == 0:
            return
        dual = _quiet(pathAlgebra.dualPathAlgebra, pathAlg)
        for vertex in allowed:
            for signed in (vertex, -vertex):
                target = pathAlg if signed > 0 else dual
                if not _quiet(mutation.mutationIsPossibleAtVertex, target, vertex):
                    continue
                nextAlg = _quiet(mutation.quiverMutationAtVertices, _copy(pathAlg), [signed])
                if nextAlg is None:
                    continue
                key = _stateKey(nextAlg)
                sequence = history + [signed]
                reached = asRelLengths(nextAlg, length)
                if reached is not None:
                    name = lines.className(reached)
                    if name != startName and (name not in best or len(sequence) < len(best[name])):
                        best[name] = sequence
                if seen.get(key, -1) >= steps - 1:
                    continue
                seen[key] = steps - 1
                walk(nextAlg, steps - 1, sequence)

    walk(nakayama.LinearNakayamaAlgebra(length, relLengths), maxSteps, [])
    return best


def smallPatterns(maxRelations = 3, maxArrows = 4, maxWidth = 7):
    """Candidate relation patterns to plant in the middle of a quiver.

    Every set of up to `maxRelations` relations, each of 2 to `maxArrows` arrows,
    with strictly increasing starts and ends -- the admissibility condition -- and
    spanning at most `maxWidth` arrows.  Normalised so the leftmost relation
    starts at 1.
    """
    import itertools

    patterns = set()
    for count in range(1, maxRelations + 1):
        for starts in itertools.combinations(range(1, maxWidth + 1), count):
            for arrows in itertools.product(range(2, maxArrows + 1), repeat = count):
                relations = list(zip(starts, arrows))
                ends = [start + arrow for start, arrow in relations]
                if any(later <= earlier for earlier, later in zip(ends, ends[1:])):
                    continue
                shift = relations[0][0] - 1
                normalised = tuple((start - shift, arrow) for start, arrow in relations)
                if patternWidth(normalised) > maxWidth:
                    continue
                patterns.add(normalised)
    return sorted(patterns)


def discoverLocalMoves(patterns, maxSteps = 3, margin = 3, embeddings = ((13, 4), (14, 5)),
                       minOccurrences = 2, progress = False):
    """Discover local rewrites by planting patterns in the middle of long quivers.

    `embeddings` is a list of (quiver length, offset) to plant each pattern at.
    Using more than one, at different lengths and offsets, is what rules out a
    rewrite that only holds because an end of the quiver happened to be nearby.

    Returns {description: [(length, before, after), ...]} for descriptions seen at
    `minOccurrences` distinct embeddings.
    """
    seen = {}
    for pattern in patterns:
        width = patternWidth(pattern)
        for length, offset in embeddings:
            relLengths = embedPattern(length, pattern, offset)
            if relLengths is None:
                continue
            if progress:
                print('  {0} in A_{1} at {2}'.format(pattern, length, offset), flush = True)
            centreLo = offset + min(s for s, _a in pattern)
            centreHi = offset + max(s + a - 1 for s, a in pattern) - 1
            reached = localMutationSequences(
                length, relLengths, centreLo, centreHi, maxSteps, margin)
            for name, sequence in reached.items():
                description = describeLink(length, relLengths, [int(c) for c in name], sequence)
                if description is None:
                    continue
                seen.setdefault(description, []).append(
                    (length, lines.className(relLengths), name))
    return {d: places for d, places in seen.items()
            if len({p[0] for p in places}) >= minOccurrences or len(places) >= minOccurrences}


# ---------------------------------------------------------------------------
# Rule families
#
# Discovery finds rules one window width at a time, so a rule that holds for
# every relation length shows up only at the widths the search happened to
# reach.  The pair slide is the clearest case: it was found at relation lengths
# 3 and 4, because those are the widths that fit at the lengths being searched,
# and it in fact holds for every length with the same two mutations.  Where a
# family is known, generate it rather than waiting for discovery to stumble on
# each member.
# ---------------------------------------------------------------------------


def pairSlideRules(maxRelationLength = 9):
    """The pair slide for every relation length, both directions.

    Two relations of equal length l starting at consecutive vertices --
    maximally overlapping -- slide one arrow along the quiver, provided no other
    relation shares an arrow with their span:

    * two *right* mutations at the first relation's source move them **left**;
    * two *left* mutations at the second relation's target move them **right**.

    The mutation count is two whatever l is; only the window widens.  Verified
    for l = 2 through 7 at four lengths each, 22 confirmations apiece with no
    failures, so the family is generated up to `maxRelationLength` rather than
    listed.
    """
    rules = []
    for relationLength in range(2, maxRelationLength + 1):
        width = relationLength + 2
        rules.append((width,
                      ((1, relationLength), (2, relationLength)),
                      ((0, relationLength), (1, relationLength)),
                      (2, 2)))
        rules.append((width,
                      ((0, relationLength), (1, relationLength)),
                      ((1, relationLength), (2, relationLength)),
                      (-width, -width)))
    return rules


def shortRelationSlideRules(maxDistance = 7):
    """A lone relation of two arrows travelling d arrows, in d mutations.

    A relation of two arrows with nothing else in its window moves d arrows
    right under the d left mutations at the window's vertices 3, 4, ..., d + 2,
    and back under the d right mutations at d + 1, d, ..., 2.  The window is
    d + 2 arrows wide.

    Unlike the pair slide, whose two mutations serve every relation length, this
    family's **mutation count grows with its parameter** -- which is why
    discovery only ever found its first three members: a search bounded at three
    mutations cannot see d >= 4, however simple the statement is.  That is
    H-008's prediction, and this is the family that confirms it (F-020).
    Verified for d = 1 to 7, both directions, at the four lengths d + 3 .. d + 6
    each -- up to A_14 and its 742900 LNAs -- with 63 confirmations apiece and no
    failures.
    """
    rules = []
    for distance in range(1, maxDistance + 1):
        width = distance + 2
        rules.append((width, ((0, 2),), ((distance, 2),),
                      tuple(-vertex for vertex in range(3, distance + 3))))
        rules.append((width, ((distance, 2),), ((0, 2),),
                      tuple(range(distance + 1, 1, -1))))
    return rules


def _withFamilies(listed):
    """The listed rules together with the generated families, deduplicated."""
    combined = list(listed)
    seen = set(combined)
    for rule in pairSlideRules() + shortRelationSlideRules():
        if rule not in seen:
            seen.add(rule)
            combined.append(rule)
    return combined


DISCOVERED_MOVES = VERIFIED_MOVES
VERIFIED_MOVES = _withFamilies(DISCOVERED_MOVES)
