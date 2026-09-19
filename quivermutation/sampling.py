"""Drawing LNAs at lengths where there are too many to enumerate.

Every classification in this repo so far has walked *all* Catalan(n - 1) LNAs of
a length.  That is 1430 rows at n = 9 and 58786 at n = 12, and 1767263190 at
n = 20, so past about n = 13 the complete pass is not the question that can be
asked.  What can be asked is a sample: draw LNAs uniformly, put each through the
cheap part of the pipeline, and count how often the answer is the interesting
one.

The point of asking at all is that **the short lengths are unrepresentative by
construction**.  In a quiver of length 8 every vertex is within three arrows of
an end, so every configuration is in reach of the anchored rules; whatever needs
a genuinely deep interior to exist cannot appear there at any depth.  A sample
at n = 16 costs the same per row as a row at n = 9 and is the only way to see
whether the leftovers thin out, hold steady, or take shapes the short lengths
have no room for.

**Uniformly, and exactly so.**  Rejection sampling over relation-length vectors
would be biased towards the sparse ones and is not what this does.  The count is
a dynamic program over the same recursion `lines.generateAllPossibleLineRelations`
enumerates by -- an LNA is a set of relations whose starts and whose ends are
both strictly increasing, each spanning at least two arrows -- and the sample is
drawn backwards through that table, so every LNA of the length has exactly the
same chance.  `countLNAs` agrees with the enumerator wherever the enumerator can
still be run, which is the test.
"""

import random


def relationSetCounts(length):
    """`counts[L][s]`: relation sets on a line of `L` vertices, last start `s`.

    `s = 0` means the set has no relations at all.  The recursion is the
    enumerator's: a relation set on `L` vertices either is one on `L - 1`
    vertices, or is one on `L - 1` vertices with a relation `i -> L` added,
    where `i` is after the previous relation's start and leaves at least two
    arrows.

    Returned as a list indexed by `L` from 0 to `length`, each entry a list
    indexed by `s` from 0 to `length`.  Small enough to hold and exact, being
    Python integers throughout.
    """
    if length < 0:
        raise ValueError("a line has a non-negative number of vertices")
    counts = [[0] * (length + 1) for _ in range(length + 1)]
    for shorter in range(0, min(2, length) + 1):
        counts[shorter][0] = 1
    for size in range(3, length + 1):
        previous = counts[size - 1]
        running = 0
        for start in range(0, size + 1):
            # `running` is the number of sets on `size - 1` vertices whose last
            # relation starts strictly before `start`, which is exactly what a
            # relation `start -> size` may be appended to.
            counts[size][start] = previous[start]
            if 1 <= start <= size - 2:
                counts[size][start] += running
            running += previous[start]
    return counts


def countLNAs(length):
    """How many LNAs there are on a line of this length: Catalan(length - 1)."""
    if length <= 0:
        return 0
    return sum(relationSetCounts(length)[length])


def sampleRelations(length, rng = None, counts = None):
    """One LNA of the length, uniformly at random, as `(start, end)` pairs.

    The pairs are in increasing order of both start and end, which is the order
    `lines` writes them in.  Pass `counts` to reuse the table across many draws;
    building it is `O(length^2)` and a draw is `O(length^2)` as well, so for a
    run of thousands it is worth passing.
    """
    if length < 3:
        return []
    rng = random.Random() if rng is None else rng
    counts = relationSetCounts(length) if counts is None else counts

    lastStart = _choose(rng, counts[length])
    relations = []
    size = length
    while size >= 3:
        previous = counts[size - 1]
        withoutHere = previous[lastStart]
        # Adding `lastStart -> size` is available only when that relation is
        # legal, and then it may be appended to any shorter set whose own last
        # relation starts strictly earlier.
        withHere = (sum(previous[:lastStart])
                    if 1 <= lastStart <= size - 2 else 0)
        if withHere and rng.randrange(withoutHere + withHere) >= withoutHere:
            relations.append((lastStart, size))
            lastStart = _choose(rng, previous[:lastStart])
        size -= 1
    relations.reverse()
    return relations


def _choose(rng, weights):
    """An index drawn with the given integer weights, exactly and without floats."""
    total = sum(weights)
    if total <= 0:
        raise ValueError("nothing to choose from")
    point = rng.randrange(total)
    for index, weight in enumerate(weights):
        point -= weight
        if point < 0:
            return index
    raise AssertionError("weights changed under the draw")


def relLengthsOf(length, relations):
    """`(start, end)` pairs -> the per-vertex relation lengths the repo uses."""
    relLengths = [0] * max(0, length - 2)
    for start, end in relations:
        relLengths[start - 1] = end - start
    return relLengths


def sampleRelLengths(length, rng = None, counts = None):
    """One LNA of the length, uniformly, as per-vertex relation lengths."""
    return relLengthsOf(length, sampleRelations(length, rng, counts))


def drawFor(length, seed, index, counts = None):
    """The `index`-th LNA of a seeded run, reproducibly and independently.

    A batch job has to be able to stop after 4000 draws and resume at the 4001st
    without keeping the generator's state, and it has to be able to re-run draw
    2317 alone to look at it again.  So each draw gets its own generator, seeded
    by the run's seed and the draw's number, rather than being the next value of
    one long stream.

    The seed is the string `"<seed>/<index>"`, not the pair: `random.Random`
    takes an int, a float, a str or bytes and refuses a tuple.
    """
    return sampleRelLengths(length, random.Random("{0}/{1}".format(seed, index)),
                            counts)


# -- what a drawn LNA is put through --------------------------------------
#
# The cheap part of the pipeline, in the order that settles the most for the
# least, so that a long run spends its time on the rows that are still open:
#
# 1. The quipu theorem names an LNA outright when its relations are almost
#    separate (F-021), and that is a string test on the relation lengths.
# 2. Otherwise the moves -- the rule table, the edge doubling, the free move and
#    the double mutation -- may carry it to one that is.  `freeMoves.orbitOf`
#    walks from this row alone rather than partitioning the length, which is the
#    only affordable direction above n = 13.
# 3. What neither settles is a **leftover**, and leftovers are the whole point
#    of sampling a long length: at n = 9, 10 and 11 they are a handful of orbits
#    and the open question about them is H-013.  Whether their rate falls, holds
#    or rises with the length, and whether they take shapes the short lengths
#    have no room for, is what a sample measures.


def probe(length, relLengths, orbitLimit = 20000):
    """Put one drawn LNA through the cheap pipeline, and say where it stopped.

    Returns a dict with the LNA's name, its overlap profile, and `settledBy`,
    one of

    * `'theorem'`  -- almost separate, so the quipu theorem names it;
    * `'moves'`    -- the moves carry it to one that is, and `movesTo` says
                      which and `orbit` how many rows were walked to find it;
    * `'leftover'` -- neither, within the orbit limit.

    A `'leftover'` is a lower bound and not a proof: `orbitOf` walks forwards
    only and stops at `orbitLimit`, so a row it does not place may still be
    placeable.  That is the same caveat the exhaustive runs carry and is why the
    limit is recorded alongside the answer.
    """
    from . import freeMoves as fm
    from . import overlap as ov

    relLengths = tuple(relLengths)
    profile = ov.overlapProfile(list(relLengths))
    record = {
        'length': length,
        'name': ''.join(str(value) for value in relLengths),
        'relations': sum(1 for value in relLengths if value),
        'maxOverlap': ov.maxOverlap(list(relLengths)),
        'overlapProfile': list(profile),
        'orbitLimit': orbitLimit,
    }
    if ov.isAlmostSeparate(length, relLengths):
        record.update(settledBy = 'theorem', orbit = 1, movesTo = record['name'])
        return record
    orbit = fm.orbitOf(length, relLengths, free = True, edges = True,
                       doubles = True, limit = orbitLimit)
    record['orbit'] = len(orbit)
    for member in sorted(orbit):
        if ov.isAlmostSeparate(length, member):
            record.update(settledBy = 'moves',
                          movesTo = ''.join(str(value) for value in member))
            return record
    record.update(settledBy = 'leftover', movesTo = None)
    return record
