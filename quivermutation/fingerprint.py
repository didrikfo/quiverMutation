"""Identifying a quiver with relations the search has met before.

`mutationSearchDepthFirst` walks mutation *sequences*, not the algebras they
reach, and the same algebra is reached along many sequences: at n = 9, depth 6,
out of the leftover orbit `3345000`, 19483 nodes are 1708 distinct algebras
(E-042).  The ratio compounds with depth -- 3.6x at depth 4, 6.3x at 5, 11.4x at
6 -- because a mutation is invertible and mutations at distant vertices commute,
so every extra level multiplies the number of routes to the same place faster
than it multiplies the places.  Everything below a repeated node is work done
twice, and this module is what lets the search not do it.

**The problem is much smaller than graph isomorphism, and it is worth saying
why.**  Vertex labels do not move under mutation (see `search.quiverKey`), so
two algebras reached from the same start are equal on the nose or not at all --
there is no relabeling to search over and nothing here is isomorphism testing.
The one thing that is genuinely ambiguous is the naming of **parallel arrows**:
the procedure hands out `networkx` edge keys in the order it builds them, so the
same algebra reached along two routes carries different keys, and NOTES.md has
recorded the want of a canonical form for that case since arrows were named.

That case is finite and, in practice, tiny.  The ambiguity is exactly a
permutation of the arrows within each bundle of parallel ones, so the canonical
form is the lexicographic least relabeling, at a cost of the product of the
bundle sizes' factorials.  Over the same two searches, every node that has
parallel arrows at all has **one** bundle of **two** arrows -- cost 2 -- and such
nodes are 0.2% to 3% of the walk (E-042).  So `canonicalKey` is exact, and the
`cap` that guards against a bundle structure nobody has yet seen has never
fired.

**The other ambiguity is the sign gauge, and it is the one that was costing
real work.**  Rescaling an arrow by -1 is an automorphism of the path algebra,
so it changes the presentation and not the algebra -- but the procedure is
sensitive to it, and mutating two presentations of one algebra at the same
vertex gives two presentations of one result.  The walk was visiting both and
walking the whole subtree below each (F-050).  Quotienting by it is linear
algebra over F_2 and costs nothing; see the note above `_gaugeGenerators`.

**What is probabilistic here is storage, not identity.**  `digest` is a 128-bit
hash of the canonical key, for a run whose visited set will not fit in memory as
keys.  It does not make identification cheaper -- the canonical key has to be
built to hash it -- it makes it *smaller*, about 30x.  Its collisions are the
usual birthday count: at 10^8 distinct algebras a 128-bit digest collides with
probability about 1.5e-23, and a 64-bit one with probability about 2.7e-4, or
one run in 3700.  Neither is a reason to use 64 bits, but the second is the
reason 64 is not absurd, and `collisionChance` is there to be asked rather than
guessed at.

**A collision costs recall and never soundness.**  A false match prunes a
subtree that was not in fact visited, so the search can *miss* a line or a
meeting point it would otherwise have found.  It cannot invent one: every
algebra reported is reported because the walk actually reached it, and a merge
certificate is a mutation path that `lnaMoves.verifyMove` re-runs.  So the
failure mode of a too-small digest is a weaker negative result, in a search
whose negative results are already only lower bounds on depth.  That asymmetry
is what makes hashing an acceptable trade here at all.

**Dedup is presentation-level, not ideal-level.**  Modulo the two ambiguities
above, the key is the presentation: two presentations that generate the same
ideal by genuinely different generators are still two keys.  That is the
conservative direction -- it costs recall, like a collision does -- and closing
it would mean a canonical basis of the ideal per source-target pair, which costs
what the Cartan matrix costs and is not obviously worth it.  See NOTES.md.
"""

import hashlib
import itertools
import math

from . import arrowPaths
from . import procedure


#: The largest number of arrow-key relabelings `canonicalKey` will minimise
#: over before giving up and calling a node un-dedupable.  The product of the
#: bundle sizes' factorials; 720 admits a bundle of six, or three of two and one
#: of three, which is far beyond anything a search has produced.
DEFAULT_CAP = 720


def arrowBundles(quiver):
    """The parallel-arrow bundles, as `(tail, head) -> sorted list of keys`.

    Only bundles of two or more arrows are returned: a bundle of one has nothing
    to permute and contributes nothing to the cost.
    """
    bundles = {}
    for tail, head, key in quiver.edges(keys = True):
        bundles.setdefault((tail, head), []).append(key)
    return {ends: sorted(keys) for ends, keys in bundles.items() if len(keys) > 1}


def relabelingCost(quiver):
    """How many arrow-key relabelings `canonicalKey` would have to try.

    The product of the bundle sizes' factorials, which is 1 -- one relabeling,
    the identity -- exactly when the quiver has no parallel arrows.
    """
    cost = 1
    for keys in arrowBundles(quiver).values():
        cost *= math.factorial(len(keys))
    return cost


def _edgeKey(quiver):
    """The quiver itself, as a sorted multiset of `(tail, head, multiplicity)`.

    Arrow keys are deliberately not in it: they are the part that is not
    canonical.  The multiplicity is, and it is what tells a parallel pair from a
    single arrow.
    """
    counts = {}
    for tail, head, _key in quiver.edges(keys = True):
        counts[(tail, head)] = counts.get((tail, head), 0) + 1
    return tuple(sorted((tail, head, count) for (tail, head), count in counts.items()))


def _normaliseSign(terms):
    """A relation and its negation are one relation, so fix which one we write.

    `terms` is a sorted tuple of `(path, coefficient)`.  The sign is chosen to
    make the first coefficient positive, which is well defined because the paths
    are sorted and a relation is never the empty combination.
    """
    if terms and terms[0][1] < 0:
        return tuple((path, -coefficient) for path, coefficient in terms)
    return terms


# -- the sign gauge -------------------------------------------------------
#
# Rescaling an arrow by a nonzero scalar is an automorphism of the path
# algebra: it changes the presentation and not the algebra.  For the signs that
# is the group {-1, +1}^arrows, acting on a term by the product of the scalars
# of the arrows the term's path runs along, and {-1, +1}^relations on top of it,
# since a relation and its negation generate the same ideal.
#
# The procedure is sensitive to that choice and the search therefore walks it.
# Research F-050: the two presentations `p = 0` and `-p = 0` of one relation are
# one algebra, and mutating them at the same vertex gives presentations that
# differ by the sign of one term -- so a walk that does not quotient by the
# gauge visits the same algebra twice, under two names, and walks the whole
# subtree below each.
#
# Quotienting is linear algebra over F_2.  Write 0 for a positive coefficient
# and 1 for a negative one, so a presentation's signs are a vector in
# F_2^terms; rescaling arrow `a` adds to it the vector that is 1 on the terms
# whose path uses `a` an odd number of times, and negating relation `r` adds
# the vector that is 1 on `r`'s terms.  The reachable presentations are exactly
# the coset `sigma + Image(M)` for the matrix `M` of those generators, so the
# canonical one is `sigma` reduced against a fixed echelon basis of Image(M).
# That is a few dozen XORs on integers and it is exact.


def _gaugeGenerators(shapes, arrows):
    """The rescalings' effect on the sign vector, as bitmasks over the terms.

    `shapes` is the relation list as `(path, magnitude)` tuples, already in the
    order the term indices follow.  One generator per arrow -- flip every term
    whose path uses it an odd number of times -- and one per relation, flipping
    all of that relation's terms.
    """
    index = {}
    for position, shape in enumerate(shapes):
        for path, _magnitude in shape:
            index[(position, path)] = len(index)

    generators = []
    for arrow in arrows:
        mask = 0
        for position, shape in enumerate(shapes):
            for path, _magnitude in shape:
                if path.count(arrow) % 2:
                    mask |= 1 << index[(position, path)]
        if mask:
            generators.append(mask)
    for position, shape in enumerate(shapes):
        mask = 0
        for path, _magnitude in shape:
            mask |= 1 << index[(position, path)]
        if mask:
            generators.append(mask)
    return generators, index


def _reduceModuloSpan(vector, generators):
    """The canonical representative of `vector + span(generators)` over F_2.

    The span is put in echelon form by leading bit, highest first, and the
    vector is reduced against it.  Two vectors in one coset reduce to the same
    thing, which is the only property wanted of it.
    """
    basis = []
    for generator in generators:
        for pivot in basis:
            generator = min(generator, generator ^ pivot)
        if generator:
            basis.append(generator)
            basis.sort(reverse = True)
    for pivot in basis:
        vector = min(vector, vector ^ pivot)
    return vector


def _relationsUnder(relations, renaming, arrows, gauge = True, cap = DEFAULT_CAP):
    """The relation set with arrow keys rewritten and the sign gauge quotiented.

    Each relation becomes a sorted tuple of `(path, coefficient)` and the
    relations are sorted, so the order they were listed in does not show.  With
    `gauge`, the signs are then replaced by the canonical representative of
    their coset under arrow rescaling -- see the note above `_gaugeGenerators`
    -- so two presentations of one algebra that differ by rescaling arrows give
    one answer.

    Relations are ordered by *shape*, meaning paths and coefficient magnitudes,
    which the gauge does not touch.  Two relations of the same shape leave that
    order ambiguous, so the answer is minimised over the ways of breaking the
    tie; `None` when there are more than `cap` of them, which no presentation a
    search has produced comes close to.  `gauge = False` is the presentation
    itself, for measuring what the quotient is worth.
    """
    written = []
    for relation in relations:
        terms = sorted(
            (tuple(renaming.get(arrow, arrow) for arrow in path), coefficient)
            for path, coefficient in relation.items()
        )
        written.append(tuple(terms))
    if not gauge:
        return tuple(sorted(_normaliseSign(terms) for terms in written))

    shapes = [tuple((path, abs(coefficient)) for path, coefficient in terms)
              for terms in written]
    groups = {}
    for position, shape in enumerate(shapes):
        groups.setdefault(shape, []).append(position)
    tied = [positions for positions in groups.values() if len(positions) > 1]
    cost = 1
    for positions in tied:
        cost *= math.factorial(len(positions))
    if cost > cap:
        return None

    fixed = sorted(groups)
    best = None
    for assignment in itertools.product(*[itertools.permutations(groups[shape])
                                          for shape in fixed]):
        order = [position for positions in assignment for position in positions]
        candidate = _gaugeCanonical([written[position] for position in order],
                                    [shapes[position] for position in order], arrows)
        if best is None or candidate < best:
            best = candidate
    return best


def _gaugeCanonical(written, shapes, arrows):
    """The relations of `written`, signed by the canonical member of their coset."""
    generators, index = _gaugeGenerators(shapes, arrows)
    vector = 0
    for position, terms in enumerate(written):
        for path, coefficient in terms:
            if coefficient < 0:
                vector |= 1 << index[(position, path)]
    vector = _reduceModuloSpan(vector, generators)
    result = []
    for position, shape in enumerate(shapes):
        result.append(tuple(
            (path, -magnitude if vector >> index[(position, path)] & 1 else magnitude)
            for path, magnitude in shape))
    return tuple(result)


def _renamings(bundles):
    """Every way of renaming the arrow keys within the bundles, as dicts.

    A bundle's arrows are renamed to `0, 1, ... ` in some order, so the result
    does not depend on which keys `networkx` happened to hand out -- only on
    which arrow of the bundle each path runs along.
    """
    if not bundles:
        yield {}
        return
    ends = sorted(bundles)
    choices = [itertools.permutations(range(len(bundles[end]))) for end in ends]
    for assignment in itertools.product(*choices):
        renaming = {}
        for end, order in zip(ends, assignment):
            tail, head = end
            for key, position in zip(bundles[end], order):
                renaming[(tail, head, key)] = (tail, head, position)
        yield renaming


def canonicalKey(pathAlg, cap = DEFAULT_CAP, gauge = True):
    """A hashable key equal for two algebras exactly when their presentations are.

    Returns `(edge multiset, relation set)`, where the relation set is written
    over arrow keys renumbered canonically within each parallel bundle -- the
    lexicographic least such writing, minimised over the bundles' permutations.
    Equal keys mean equal algebras, labels and all; different keys mean the
    presentations differ, which is *almost* the same as the algebras differing
    and is conservative where it is not (see the module docstring).

    `None` when the quiver would need more than `cap` relabelings, which no
    search has yet produced.  A caller that gets `None` must treat the node as
    one it has not seen before -- that is the safe direction -- which is what
    `Visited.seen` does.

    Where the quiver has no parallel arrows there is one relabeling and this is
    `search.quiverKey` with the coefficients added and the sign gauge taken out:
    exact, and no more expensive than reading the relations off.

    `gauge = False` keys the presentation as written, without quotienting by
    arrow rescaling.  That is a strictly finer partition and is for measuring
    what the quotient is worth, not for a search to walk on.
    """
    bundles = arrowBundles(pathAlg.quiver)
    if bundles:
        cost = 1
        for keys in bundles.values():
            cost *= math.factorial(len(keys))
        if cost > cap:
            return None
    relations = procedure.relationsFrom(pathAlg)
    best = None
    for renaming in _renamings(bundles):
        arrows = [renaming.get(arrow, arrow)
                  for arrow in arrowPaths.arrowsOf(pathAlg.quiver)]
        candidate = _relationsUnder(relations, renaming, arrows,
                                    gauge = gauge, cap = cap)
        if candidate is None:
            return None
        if best is None or candidate < best:
            best = candidate
    return (_edgeKey(pathAlg.quiver), best)


def digest(key, digestBits = 128):
    """A `digestBits`-wide hash of a canonical key, as an int.

    For a visited set too large to hold as keys.  `None` in, `None` out, so a
    node `canonicalKey` refused stays refused.  See the module docstring for the
    collision arithmetic and for why a collision costs only recall.
    """
    if key is None:
        return None
    if digestBits % 8 or not 32 <= digestBits <= 512:
        raise ValueError("digestBits must be a multiple of 8 between 32 and 512")
    raw = hashlib.blake2b(repr(key).encode("utf-8"), digest_size = digestBits // 8)
    return int.from_bytes(raw.digest(), "big")


def collisionChance(count, digestBits = 128):
    """The chance that `count` distinct keys collide somewhere, approximately.

    The birthday bound `1 - exp(-count^2 / 2^(bits+1))`, for choosing a width
    against the size of a run rather than by superstition.
    """
    if count < 2:
        return 0.0
    return -math.expm1(-(count * (count - 1)) / 2.0 ** (digestBits + 1))


class Visited:
    """The set of algebras a search has already expanded, and at what depth.

    A plain visited set is wrong for a depth-bounded walk: a node first met with
    one mutation left and met again with four still has three levels below it
    that were never walked.  So what is stored is the largest *remaining* depth
    the node has been expanded at, and `seen` says to skip only a node that has
    already been expanded at least as deep.  That makes the deduped walk reach
    exactly what the plain walk reaches (E-043 checks this, node for node, over
    every LNA of n = 5 to 8).

    What it does *not* preserve is which mutation path is recorded first: the
    plain walk reaches a node once per route and the deduped walk once, so a
    caller keeping the shortest path may be handed a longer one.  Nothing in the
    repo depends on the path being shortest -- `quiversReachedFrom` keeps the
    shortest it is *shown*, and a merge certificate is checked, not minimised --
    but a caller that did would have to search breadth-first instead.

    With `digestBits` the keys are stored as hashes, about 30x smaller, at the
    cost described in the module docstring.  With `cap` raised or lowered, the
    point at which a parallel-arrow node is called un-dedupable moves; an
    un-dedupable node is always walked, never skipped.
    """

    def __init__(self, digestBits = None, cap = DEFAULT_CAP, gauge = True):
        self.digestBits = digestBits
        self.cap = cap
        self.gauge = gauge
        self.depths = {}
        self.hits = 0
        self.misses = 0
        self.refused = 0

    def _key(self, pathAlg):
        key = canonicalKey(pathAlg, cap = self.cap, gauge = self.gauge)
        if key is None:
            return None
        return key if self.digestBits is None else digest(key, self.digestBits)

    def seen(self, pathAlg, remainingDepth):
        """Whether this algebra has already been expanded at least this deep.

        Records the visit when the answer is no, so a caller asks once per node
        and does not have to record separately.
        """
        key = self._key(pathAlg)
        if key is None:
            self.refused += 1
            return False
        previous = self.depths.get(key)
        if previous is not None and previous >= remainingDepth:
            self.hits += 1
            return True
        self.depths[key] = remainingDepth
        self.misses += 1
        return False

    def summarise(self):
        """What the dedup did, for a run to print and a checkpoint to record."""
        visits = self.hits + self.misses + self.refused
        return {
            'nodes': visits,
            'distinct': len(self.depths),
            'skipped': self.hits,
            'refused': self.refused,
            'ratio': (visits / len(self.depths)) if self.depths else 1.0,
            'digestBits': self.digestBits,
        }

    def __len__(self):
        return len(self.depths)
