"""The doubling at an end of the quiver, which the window encoding cannot state.

A rule in `lnaMoves` is a rewrite on a window of arrows, and `matchesAt` insists
that every relation meeting the window lies entirely inside it -- otherwise the
rewrite would depend on something it does not describe.  That is the right
default and it is what keeps the table honest (R-009).  It also makes one true
family inexpressible, and research F-029 is that family:

    A relation of `l` arrows at the source, with nothing starting at the
    vertices 2 .. l, gains a second relation of `l` arrows at vertex 2 under two
    left mutations at vertex `l + 1`.

Its window is the arrows 1 .. l + 1, and a relation of the quiver is allowed to
*start* on the window's last arrow and run out the far side -- the rewrite never
touches it.  `describeLink` refuses to describe that, so no discovery run could
ever have listed this rule, which is why 1794 verified rules do not contain it.

The dual is the same thing at the sink, by the relation dual of F-026: a
relation of `l` arrows ending at the sink, with nothing ending at the `l - 1`
vertices before it, gains one ending at vertex `n - 1` under two right mutations
at vertex `n - l`.

Both directions of both are checked against the mutation engine in
`tests/test_edge_moves.py`, the same three ways `lnaMoves.verifyMove` checks a
rule: the predicted LNA, every mutation admissible, and the Coxeter polynomial
held fixed.
"""

from . import lnaMoves


def _admissible(length, relLengths):
    """Direct, not a lookup: membership of `allRelationLengths` would rebuild a
    set of every LNA of the length on each call, which at n = 12 is 58786 rows
    per rewrite."""
    return lnaMoves.isAdmissible(length, list(relLengths))


def sourceDoubling(length, relLengths):
    """(doubled relation lengths, mutation sequence), or None.

    The condition is that a relation of `l` arrows starts at vertex 1 and
    nothing starts at the vertices 2 .. l.  A relation starting at vertex
    `l + 1` -- on the window's last arrow -- is allowed, and is what the window
    encoding cannot express.
    """
    relLengths = list(relLengths)
    arrows = relLengths[0] if relLengths else 0
    if arrows < 2 or len(relLengths) < 2:
        return None
    if any(relLengths[1:arrows]):
        return None
    result = list(relLengths)
    result[1] = arrows
    if not _admissible(length, result):
        return None
    return tuple(result), [-(arrows + 1), -(arrows + 1)]


def sourceCollapse(length, relLengths):
    """The inverse: two equal relations at vertices 1 and 2 become one.

    Defined as the exact inverse of `sourceDoubling` -- the LNA whose doubling
    is this one -- rather than by inverting its condition by hand, which is how
    the first version of `sinkCollapse` came to fire on only a third of the
    cases it should have.
    """
    relLengths = tuple(relLengths)
    if len(relLengths) < 2 or relLengths[1] == 0:
        return None
    candidate = list(relLengths)
    candidate[1] = 0
    doubled = sourceDoubling(length, candidate)
    if doubled is None or doubled[0] != relLengths:
        return None
    # Two right mutations at the *source*, not at the vertex the doubling
    # mutated: the procedure relabels, so a left mutation is not undone by a
    # right mutation at the same vertex.  This is the sequence
    # `endPairCollapseRules` already uses for the case with no spectator.
    return tuple(candidate), [1, 1]


def sinkDoubling(length, relLengths):
    """The dual of `sourceDoubling`, against the sink."""
    relLengths = list(relLengths)
    position = None
    for start, arrows in lnaMoves.relationsOf(relLengths):
        if start + arrows == length:
            position, count = start - 1, arrows
            break
    if position is None or count < 2:
        return None
    ends = [start + arrows for start, arrows in lnaMoves.relationsOf(relLengths)]
    if any(length - count + 1 <= end <= length - 1 for end in ends):
        return None
    target = position - 1
    if target < 0 or relLengths[target]:
        return None
    result = list(relLengths)
    result[target] = count
    if not _admissible(length, result):
        return None
    return tuple(result), [length - count, length - count]


def sinkCollapse(length, relLengths):
    """The inverse of `sinkDoubling`, defined the same way as `sourceCollapse`."""
    relLengths = tuple(relLengths)
    for start, arrows in lnaMoves.relationsOf(list(relLengths)):
        if start + arrows != length - 1:
            continue
        candidate = list(relLengths)
        candidate[start - 1] = 0
        doubled = sinkDoubling(length, candidate)
        if doubled is not None and doubled[0] == relLengths:
            # The dual of the source collapse: two left mutations at the sink.
            return tuple(candidate), [-length, -length]
    return None


MOVES = (sourceDoubling, sourceCollapse, sinkDoubling, sinkCollapse)


def rewritesOf(length, relLengths):
    """Every edge move out of an LNA, as (relation lengths, sequence) pairs."""
    reached = []
    for move in MOVES:
        answer = move(length, list(relLengths))
        if answer is not None:
            reached.append(answer)
    return reached
