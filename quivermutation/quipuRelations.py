"""Quipu quivers that do carry relations, against the LNAs of the same length.

The quipu theorem is a statement about quipu quivers with **no** relations: an
LNA with almost separate relations is derived equivalent to the path algebra of
a quipu, and that algebra is hereditary.  Nothing says the quipu shape stops
being special once relations are put on it, and relations are what the mutation
procedure spends its time on -- a walk from an LNA to its quipu passes through
quivers that are neither.  So the question here is the other half of the one
`treeSearch` asks:

    is some quipu quiver *with relations* derived equivalent to a linear
    Nakayama algebra that lies in no quipu class?

The enumeration.  For each quipu of the order, each orientation of its edges
(up to the tree's automorphisms, which permute orientations without changing the
quiver up to isomorphism), and each admissible monomial ideal on the result.  A
tree has no oriented cycles, so *any* set of directed paths of two or more
arrows generates an admissible ideal, and the ideal is determined by its minimal
generators -- an antichain under "is a contiguous subpath of".  Antichains are
what `relationSets` walks.

`minArrows` bounds the shortest relation allowed.  With `minArrows = 2` the
enumeration is everything; with 3 it leaves out every ideal that has a relation
of two arrows, which `corollary:lengthtworelations` of arXiv:2310.08346 says is
free *for a Nakayama algebra*.  That corollary is not stated for quipus, so the
restriction is an assumption and not a theorem here -- `freeRelationCheck`
measures how far it holds on this family, and the counts it gives are the reason
`minArrows = 3` is the default at the orders where everything is too much.

The cost.  At order 11 there are 7.2 million ideals with every relation of three
arrows or more, which rules out computing a Coxeter polynomial per algebra the
symbolic way.  Two things make it affordable.  The Cartan matrix of a monomial
algebra on a tree is a **bitmask**: bit `u * n + v` says the path from u to v
survives, and walking the antichain tree maintains it in one OR per step.  And
the polynomial only has to be computed for the algebras that could match, which
a single number decides: `det(2 C^T + C)` is the polynomial at 2, so an algebra
whose value is not one of the values the LNAs of the length take cannot match
anything, and that determinant is a batched float one over thousands of algebras
at a time.  Everything that survives the filter is recomputed exactly, in
integers, so the float is a sieve and never an answer.
"""

import collections
import functools
import itertools

import networkx as nx
import numpy as np

from . import coxeterTables
from . import invariants
from . import pathAlgebra
from . import quipuForms


FINGERPRINT_POINT = 2
BATCH = 20000


def orientedQuipus(order):
    """Every quipu of the order, oriented every way, up to its automorphisms.

    Yields `(parameters, edges, orientation, automorphisms)`, where `edges` is
    the tree's edge list on vertices `0 .. order - 1`, `orientation` a tuple of
    bits -- bit i meaning edge i runs from its first vertex to its second -- and
    `automorphisms` the underlying tree's automorphism group, which is what
    decides when two algebras on this quiver are the same one.

    Orientations in one automorphism orbit give isomorphic quivers, so only the
    lexicographically smallest of each orbit is yielded.  For a quipu with a
    symmetry -- the main string reversing onto itself, two cords of equal length
    at the same foot -- that halves the work or better; for one with none it
    yields all 2^(n-1).
    """
    for parameters in quipuForms.allQuipusOfOrder(order):
        graph = nx.convert_node_labels_to_integers(
            quipuForms.graphFromQuipuParameters(*parameters), ordering = "sorted")
        edges = sorted(tuple(sorted(edge)) for edge in graph.edges)
        edgeIndex = {edge: index for index, edge in enumerate(edges)}
        automorphisms = _automorphisms(graph)
        seen = set()
        for orientation in itertools.product((0, 1), repeat = len(edges)):
            if orientation in seen:
                continue
            orbit = {_mapOrientation(orientation, edges, edgeIndex, mapping)
                     for mapping in automorphisms}
            seen.update(orbit)
            if orientation != min(orbit):
                continue
            yield parameters, edges, orientation, automorphisms


def _automorphisms(graph):
    """The automorphism group of a small tree, as a list of vertex maps."""
    matcher = nx.algorithms.isomorphism.GraphMatcher(graph, graph)
    return [dict(mapping) for mapping in matcher.isomorphisms_iter()]


def _mapOrientation(orientation, edges, edgeIndex, mapping):
    """The orientation carried across by an automorphism of the tree."""
    mapped = [0] * len(edges)
    for index, (first, second) in enumerate(edges):
        image = (mapping[first], mapping[second])
        target = tuple(sorted(image))
        mapped[edgeIndex[target]] = (
            orientation[index] if image == target else 1 - orientation[index])
    return tuple(mapped)


def successors(order, edges, orientation):
    """The directed adjacency of one oriented tree, as lists."""
    succ = [[] for _ in range(order)]
    for (first, second), bit in zip(edges, orientation):
        if bit:
            succ[first].append(second)
        else:
            succ[second].append(first)
    return succ


def directedPaths(order, succ):
    """Every directed path of one arrow or more, as tuples of vertices."""
    found = []
    stack = [(vertex,) for vertex in range(order)]
    while stack:
        path = stack.pop()
        for onward in succ[path[-1]]:
            extended = path + (onward,)
            found.append(extended)
            stack.append(extended)
    return found


def _contains(outer, inner):
    """Whether `inner` is a contiguous subpath of `outer`."""
    span = len(inner)
    return any(outer[start:start + span] == inner
               for start in range(len(outer) - span + 1))


def relationData(order, paths, minArrows):
    """What the antichain walk needs, for one oriented quiver.

    Returns `(baseMask, candidates, killMasks, comparable)`:

    * `baseMask` -- the Cartan matrix of the quiver with no relations, as a bit
      per ordered pair, bit `u * order + v` set when the path u -> v exists.
      The diagonal is set, since e_i A e_i is always one-dimensional here.
    * `candidates` -- the paths long enough to be a relation.
    * `killMasks[i]` -- the pairs whose path contains candidate i, so choosing it
      as a relation clears exactly those bits.
    * `comparable[i]` -- a bitmask over candidate indices, those that contain
      candidate i or are contained in it, which an antichain may not use with it.
    """
    baseMask = 0
    for vertex in range(order):
        baseMask |= 1 << (vertex * order + vertex)
    for path in paths:
        baseMask |= 1 << (path[0] * order + path[-1])
    candidates = [path for path in paths if len(path) - 1 >= minArrows]
    killMasks = []
    for candidate in candidates:
        mask = 0
        for path in paths:
            if _contains(path, candidate):
                mask |= 1 << (path[0] * order + path[-1])
        killMasks.append(mask)
    comparable = []
    for index, candidate in enumerate(candidates):
        mask = 0
        for other, otherPath in enumerate(candidates):
            if other != index and (_contains(candidate, otherPath)
                                   or _contains(otherPath, candidate)):
                mask |= 1 << other
        comparable.append(mask)
    return baseMask, candidates, killMasks, comparable


def relationSets(baseMask, killMasks, comparable):
    """Every admissible ideal, as `(cartan mask, chosen candidate indices)`.

    The walk is over antichains of the containment order, taken in index order,
    so each ideal is produced exactly once, the empty one included -- that is the
    hereditary quipu itself, which is the quipu theorem's case and is kept so the
    enumeration can be checked against it.
    """
    count = len(killMasks)

    def walk(start, allowed, killed, chosen):
        yield baseMask & ~killed, chosen
        for index in range(start, count):
            if allowed >> index & 1:
                yield from walk(index + 1,
                                allowed & ~comparable[index],
                                killed | killMasks[index],
                                chosen + (index,))

    yield from walk(0, (1 << count) - 1, 0, ())


def cartanFromMask(mask, order):
    """The Cartan matrix a bitmask stands for, as nested lists of ints."""
    return [[(mask >> (source * order + target)) & 1 for source in range(order)]
            for target in range(order)]


def _fingerprints(masks, order, point = FINGERPRINT_POINT):
    """The polynomial at `point` for a batch of Cartan masks, as int64.

    `det(point * C^T + C)` is the Coxeter polynomial evaluated at `point`, by the
    identity in `invariants`.  Computed in float64 over the whole batch at once
    and rounded: the entries are at most `point + 1` and the order at most a
    dozen, so Hadamard bounds the determinant well inside float64's exact range
    and the rounding is safe.  It is a sieve regardless -- everything it selects
    is recomputed in exact integer arithmetic.
    """
    width = (order * order + 7) // 8
    buffer = b"".join(mask.to_bytes(width, "little") for mask in masks)
    bits = np.unpackbits(
        np.frombuffer(buffer, dtype = np.uint8).reshape(len(masks), width),
        axis = 1, bitorder = "little")
    transposed = bits[:, :order * order].reshape(-1, order, order).astype(np.float64)
    matrices = point * transposed + transposed.transpose(0, 2, 1)
    return np.rint(np.linalg.det(matrices)).astype(np.int64)


def _bitMatrices(masks, order):
    """A batch of Cartan masks as an (algebras, order, order) array of C^T."""
    width = (order * order + 7) // 8
    buffer = b"".join(mask.to_bytes(width, "little") for mask in masks)
    bits = np.unpackbits(
        np.frombuffer(buffer, dtype = np.uint8).reshape(len(masks), width),
        axis = 1, bitorder = "little")
    return bits[:, :order * order].reshape(-1, order, order).astype(np.int64)


# Three primes just under 2**20.  Their product is 1.15e18, and the Coxeter
# polynomial of one of these algebras at a point 0 <= t <= order is a determinant
# of a matrix with entries at most order + 1, which Hadamard bounds by
# ((order + 1) * sqrt(order))**order -- 4.0e17 at order 11.  Two such values
# therefore differ by less than the product of the primes, so *equal residues
# modulo all three means equal integers*, and the comparison needs no
# reconstruction.  The primes are small enough that a product of two residues
# stays inside int64 with room to spare, which is what lets the elimination run
# over a whole batch of matrices at once in numpy.
MODULI = (1048573, 1048571, 1048559)


def modulusBound(order):
    """Hadamard's bound on the determinants, for checking the moduli suffice."""
    return ((order + 1) * order**0.5)**order


@functools.lru_cache(maxsize = None)
def _inverseTable(prime):
    """Modular inverses of every residue, as an array to index with.

    A table rather than a `pow` per pivot: the elimination needs one inverse per
    matrix per step, and looking them up is what keeps the step vectorised.
    Entry 0 is 0, which is never used -- a zero pivot is handled before it.
    """
    table = [0] * prime
    if prime > 1:
        table[1] = 1
    for value in range(2, prime):
        table[value] = (-(prime // value) * table[prime % value]) % prime
    return np.array(table, dtype = np.int64)


def exactValues(masks, order):
    """The Coxeter polynomial of each algebra at 0, 1, ..., order, as residues.

    Two polynomials of degree at most `order` are equal exactly when they agree
    at `order + 1` points, so these values decide a match without interpolation
    and without symbols.  They are held modulo `MODULI` rather than as integers,
    for the reason given there -- equal residues mean equal values here -- and
    every step is integer arithmetic over the whole batch at once.

    Returns an array of shape (algebras, (order + 1) * len(MODULI)).  Compare it
    with what `targetValueVectors` produces and nothing else: the residues are
    not the values.
    """
    if 2 * modulusBound(order) >= _moduliProduct():
        raise ValueError("order {0} needs more moduli than {1} to be decided "
                         "this way".format(order, len(MODULI)))
    transposed = _bitMatrices(masks, order)
    cartan = transposed.transpose(0, 2, 1)
    columns = []
    for point in range(order + 1):
        matrices = point * transposed + cartan
        for prime in MODULI:
            columns.append(_modularDeterminants(matrices, prime))
    return np.stack(columns, axis = 1)


def _moduliProduct():
    product = 1
    for prime in MODULI:
        product *= prime
    return product


def _modularDeterminants(matrices, prime):
    """The determinant of every matrix in a batch, modulo one prime.

    Ordinary elimination, not Bareiss: fraction-free elimination keeps integers
    but squares their size at every step, which is exactly what overflows here,
    while modulo a prime the pivot can simply be inverted and everything stays
    below the prime.  A zero pivot needs a row swap in that one matrix, which is
    the only part that is not done to the whole batch at once.
    """
    working = matrices % prime
    count, size, _ = working.shape
    inverses = _inverseTable(prime)
    determinants = np.ones(count, dtype = np.int64)
    for step in range(size):
        pivots = working[:, step, step]
        for index in np.nonzero(pivots == 0)[0]:
            below = working[index, step + 1:, step]
            nonzero = np.nonzero(below)[0]
            if not len(nonzero):
                determinants[index] = 0
                continue
            other = step + 1 + nonzero[0]
            working[index, [step, other]] = working[index, [other, step]]
            determinants[index] = (-determinants[index]) % prime
        pivots = working[:, step, step]
        determinants = (determinants * pivots) % prime
        if step == size - 1:
            break
        scaled = inverses[pivots]
        working[:, step, step:] = (working[:, step, step:] * scaled[:, None]) % prime
        factors = working[:, step + 1:, step].copy()
        working[:, step + 1:, step:] = (
            working[:, step + 1:, step:]
            - factors[:, :, None] * working[:, step, None, step:]) % prime
    return determinants


def targetValueVectors(length, statuses = None):
    """The LNA polynomials to match against, as `{values at 0..n: coefficients}`.

    The same table `targetFingerprints` sieves with, in the form `exactValues`
    produces, so a match is an exact equality of polynomials and not a test that
    two numbers agree.
    """
    status = coxeterTables.lnaStatus(length)
    wanted = {}
    for relLengths, key in coxeterTables.lnaKeys(length).items():
        if statuses is not None and status[relLengths] not in statuses:
            continue
        values = tuple(int(sum(coefficient * point**degree
                               for degree, coefficient in enumerate(key))) % prime
                       for point in range(length + 1) for prime in MODULI)
        wanted[values] = key
    return wanted


def targetFingerprints(length, statuses = None):
    """The LNA polynomials to sieve against, as `{value at the point: keys}`.

    With `statuses` given, only the LNAs of those statuses are included, which is
    how a run asks for "LNAs in no quipu class" and nothing else.
    """
    status = coxeterTables.lnaStatus(length)
    wanted = {}
    for relLengths, key in coxeterTables.lnaKeys(length).items():
        if statuses is not None and status[relLengths] not in statuses:
            continue
        value = int(sum(coefficient * FINGERPRINT_POINT**degree
                        for degree, coefficient in enumerate(key)))
        wanted.setdefault(value, set()).add(key)
    return wanted


def isLinearlyOriented(order, edges, orientation):
    """Whether the quiver is the line 1 -> 2 -> ... -> n, i.e. an LNA's quiver.

    A quipu with no cords is the line `A_n`, and one of its orientations is the
    linear one, so the enumeration contains every LNA of the length as a member.
    Those match themselves, which is not a finding; every *other* orientation of
    the line is a genuinely different algebra and stays in.
    """
    succ = successors(order, edges, orientation)
    if any(len(onward) > 1 for onward in succ):
        return False
    hasPredecessor = [False] * order
    for vertex in range(order):
        for onward in succ[vertex]:
            hasPredecessor[onward] = True
    sources = [vertex for vertex in range(order) if not hasPredecessor[vertex]]
    if len(sources) != 1:
        return False
    visited = 1
    current = sources[0]
    while succ[current]:
        current = succ[current][0]
        visited += 1
    return visited == order


def canonicalAlgebra(order, edges, orientation, relations, automorphisms):
    """A canonical form of the quiver with relations, up to tree automorphism.

    Two ideals on the same oriented quipu can be carried onto each other by an
    automorphism of the underlying tree, in which case the algebras are
    isomorphic and one of them is redundant.  The certificate is the smallest
    (arrows, relations) pair over the group, arrows and relations both sorted.
    """
    arrows = []
    for (first, second), bit in zip(edges, orientation):
        arrows.append((first, second) if bit else (second, first))
    best = None
    for mapping in automorphisms:
        image = (tuple(sorted((mapping[tail], mapping[head]) for tail, head in arrows)),
                 tuple(sorted(tuple(mapping[vertex] for vertex in path)
                              for path in relations)))
        if best is None or image < best:
            best = image
    return best


def search(order, minArrows = 3, statuses = None, includeHereditary = False,
           includeLines = False, keepPerKey = 20, dedupe = False, progress = None):
    """Walk every quipu with relations of the order and report what matches.

    `statuses` restricts which LNAs count as a match -- pass
    `(coxeterTables.NOT_QUIPU, coxeterTables.UNPLACED)` for the question this
    module is about, or leave it None to match against every LNA of the length.
    `includeHereditary` keeps the empty ideal, which is the quipu theorem's own
    case; `includeLines` keeps the algebras whose quiver is the linearly oriented
    line, which *are* the LNAs and so match themselves.  Both are off by default
    because both are known answers rather than findings.

    `keepPerKey` bounds how many matching algebras are kept per Coxeter
    polynomial -- the simplest ones, fewest and shortest relations first -- and
    everything else is counted and dropped.  That bound is what makes the run
    affordable in memory: at order 10 the matches run to three hundred thousand
    and at order 11 to millions, and holding them all as records costs gigabytes
    to say something a count says better.  Pass None to keep every one.

    `dedupe` counts each algebra once up to isomorphism rather than once per
    orientation the enumeration reached it by, at the cost of holding every
    match's certificate.  It is what a run that wants to intersect the matches
    with a mutation class needs, and `certificates` in the result is that set.

    Returns a dict with the counts of what was walked, `counts` and `shapes` per
    polynomial, and `examples` -- the kept matches, also flattened into
    `matches`.  Each match carries the quipu, the orientation, the relations as
    vertex paths, the exact polynomial and the LNAs that share it.
    """
    targets = targetFingerprints(order, statuses)
    index = coxeterTables.lnaKeyIndex(order)
    status = coxeterTables.lnaStatus(order)
    state = {
        'order': order,
        'minArrows': minArrows,
        'keepPerKey': keepPerKey,
        'walked': 0,
        'sieved': 0,
        'matched': 0,
        'counts': collections.Counter(),
        'shapes': collections.defaultdict(collections.Counter),
        'examples': collections.defaultdict(list),
        'seen': set() if dedupe else None,
        'targetArray': np.array(sorted(targets), dtype = np.int64),
        'vectors': targetValueVectors(order, statuses),
        'index': index,
        'status': status,
    }
    for parameters, edges, orientation, automorphisms in orientedQuipus(order):
        if not includeLines and isLinearlyOriented(order, edges, orientation):
            continue
        succ = successors(order, edges, orientation)
        paths = directedPaths(order, succ)
        baseMask, candidates, killMasks, comparable = relationData(order, paths, minArrows)
        batch = []
        for mask, chosen in relationSets(baseMask, killMasks, comparable):
            if not chosen and not includeHereditary:
                continue
            state['walked'] += 1
            batch.append((mask, chosen))
            if len(batch) >= BATCH:
                _drain(batch, state, parameters, edges, orientation, automorphisms, candidates)
                batch = []
        if batch:
            _drain(batch, state, parameters, edges, orientation, automorphisms, candidates)
        if progress is not None:
            progress(parameters, orientation, state['walked'], state['matched'])
    matches = [match for group in state['examples'].values() for match in group]
    return {
        'certificates': state['seen'],
        'order': order,
        'minArrows': minArrows,
        'walked': state['walked'],
        'sieved': state['sieved'],
        'matched': state['matched'],
        'counts': dict(state['counts']),
        'shapes': {key: dict(names) for key, names in state['shapes'].items()},
        'examples': {key: list(group) for key, group in state['examples'].items()},
        'matches': matches,
    }


def _confirm(masks, state):
    """The polynomial of each sieved algebra, or None where it is not a target.

    Batched and exact while `MODULI` decides the order -- up to 11 as they stand
    -- and beyond that each one goes through `invariants` one at a time, which is
    slower and has no bound to respect.
    """
    order = state['order']
    if 2 * modulusBound(order) < _moduliProduct():
        values = exactValues(masks, order)
        return [state['vectors'].get(tuple(int(value) for value in row)) for row in values]
    keys = set(state['vectors'].values())
    found = []
    for mask in masks:
        key = invariants.coxeterCoefficients(cartanFromMask(mask, order))
        found.append(key if key in keys else None)
    return found


def _complexity(match):
    """How simple a match is, for keeping the simplest ones: relations, then arrows."""
    return (len(match['relations']), sum(len(path) for path in match['relations']))


def _drain(batch, state, parameters, edges, orientation, automorphisms, candidates):
    """Sieve one batch, count the exact matches and keep the simplest of them.

    Two stages, and the first one is only an optimisation: the float determinant
    at 2 throws out the algebras whose polynomial cannot be a target, and
    everything it keeps is then decided exactly, by the integer values of the
    polynomial at `order + 1` points.  Nothing is recorded on the strength of a
    float.
    """
    order = state['order']
    values = _fingerprints([mask for mask, _ in batch], order)
    survivors = np.nonzero(np.isin(values, state['targetArray']))[0]
    if not len(survivors):
        return
    state['sieved'] += len(survivors)
    confirmed = _confirm([batch[int(position)][0] for position in survivors], state)
    for row, position in enumerate(survivors):
        key = confirmed[row]
        if key is None:
            continue
        mask, chosen = batch[int(position)]
        relations = tuple(candidates[i] for i in chosen)
        certificate = None
        if state['seen'] is not None:
            certificate = canonicalAlgebra(order, edges, orientation, relations, automorphisms)
            if certificate in state['seen']:
                continue
            state['seen'].add(certificate)
        state['matched'] += 1
        state['counts'][key] += 1
        state['shapes'][key][quipuForms.formatQuipu(parameters)] += 1
        kept = state['examples'][key]
        keepPerKey = state['keepPerKey']
        if keepPerKey == 0:
            continue
        if keepPerKey is not None and len(kept) >= keepPerKey:
            worst = max(range(len(kept)), key = lambda i: _complexity(kept[i]))
            candidate = (len(relations), sum(len(path) for path in relations))
            if candidate >= _complexity(kept[worst]):
                continue
            kept.pop(worst)
        match = {
            'quipu': quipuForms.formatQuipu(parameters),
            'parameters': parameters,
            'edges': tuple(edges),
            'orientation': tuple(orientation),
            'relations': relations,
            'key': key,
            'lnas': tuple((coxeterTables.className(relLengths), state['status'][relLengths])
                          for relLengths in state['index'][key]),
        }
        if certificate is not None:
            match['certificate'] = certificate
        kept.append(match)


def algebraFromMatch(match):
    """The path algebra a match names, so it can be mutated or drawn.

    Vertices are numbered from 1, as everywhere else in the repo; the relations
    are the zero relations on the paths the match records.
    """
    algebra = pathAlgebra.PathAlgebra()
    order = max(max(edge) for edge in match['edges']) + 1
    algebra.add_vertices_from(range(1, order + 1))
    for (first, second), bit in zip(match['edges'], match['orientation']):
        if bit:
            algebra.add_arrow(first + 1, second + 1)
        else:
            algebra.add_arrow(second + 1, first + 1)
    for path in match['relations']:
        algebra.add_rel([[vertex + 1 for vertex in path]])
    return algebra


def freeRelationCheck(order, minArrows = 2, limit = None):
    """Does deleting a relation of two arrows keep the Coxeter polynomial here?

    `corollary:lengthtworelations` says it does for a Nakayama algebra, and
    `minArrows = 3` leans on the same thing holding for quipus.  This tests the
    consequence that is cheap to test: walk the ideals that have a relation of
    exactly two arrows, delete every such relation, and compare the polynomial
    before and after.  Returns `(checked, agreeing, disagreeing examples)`.

    A disagreement is decisive -- the polynomial is a derived invariant, so two
    algebras with different ones are not derived equivalent -- and would say the
    corollary does not extend to this family.  Agreement is evidence and no
    more.
    """
    checked = 0
    agreeing = 0
    failures = []
    for parameters, edges, orientation, _automorphisms in orientedQuipus(order):
        succ = successors(order, edges, orientation)
        paths = directedPaths(order, succ)
        baseMask, candidates, killMasks, comparable = relationData(order, paths, minArrows)
        for mask, chosen in relationSets(baseMask, killMasks, comparable):
            shortOnes = [i for i in chosen if len(candidates[i]) - 1 == 2]
            if not shortOnes:
                continue
            checked += 1
            kept = [i for i in chosen if i not in shortOnes]
            strippedMask = baseMask
            for i in kept:
                strippedMask &= ~killMasks[i]
            before = invariants.coxeterCoefficients(cartanFromMask(mask, order))
            after = invariants.coxeterCoefficients(cartanFromMask(strippedMask, order))
            if before == after:
                agreeing += 1
            else:
                failures.append({
                    'quipu': quipuForms.formatQuipu(parameters),
                    'orientation': tuple(orientation),
                    'relations': tuple(candidates[i] for i in chosen),
                    'before': before,
                    'after': after,
                })
            if limit is not None and checked >= limit:
                return checked, agreeing, failures
    return checked, agreeing, failures


# -- from the other side: which of these a class actually passes through ---
#
# A shared Coxeter polynomial is a necessary condition and no more, so the
# enumeration above produces leads, not results.  What settles one is a mutation
# path, and the cheap way to get many at once is to walk the mutation graph out
# of the LNA rather than out of each lead: every quiver a walk reaches is in the
# class by construction, so collecting the ones whose quiver is a quipu gives a
# set of *confirmed* members to intersect with the leads.
#
# The certificate has to be the same one `search` deduplicates by, which means
# canonicalising a quiver that arrives with whatever vertex numbering a mutation
# left it with.  `certificate` does that by matching the quiver's underlying tree
# against the reference quipu of its parameters -- the same graph `orientedQuipus`
# enumerates orientations of -- and taking the smallest image over every
# isomorphism between them.


def certificate(pathAlg):
    """The canonical form of a quipu algebra with monomial relations, or None.

    None when the quiver is not a quipu, when it has parallel arrows, or when any
    relation is a sum of paths rather than a single one -- all three happen along
    a mutation walk, and none of them is in the family `search` enumerates, so
    None is "outside the family" rather than a failure.
    """
    graph = quipuForms.underlyingGraph(pathAlg)
    if pathAlg.quiver.number_of_edges() != graph.number_of_edges():
        return None
    parameters = quipuForms.quipuParameters(graph)
    if parameters is None:
        return None
    if any(len(rel) != 1 for rel in pathAlg.rels):
        return None
    reference = nx.convert_node_labels_to_integers(
        quipuForms.graphFromQuipuParameters(*parameters), ordering = "sorted")
    matcher = nx.algorithms.isomorphism.GraphMatcher(graph, reference)
    best = None
    for mapping in matcher.isomorphisms_iter():
        image = (tuple(sorted((mapping[tail], mapping[head])
                              for tail, head in pathAlg.quiver.edges())),
                 tuple(sorted(tuple(mapping[vertex] for vertex in rel[0])
                              for rel in pathAlg.rels)))
        if best is None or image < best:
            best = image
    return best


def reachedQuipuAlgebras(pathAlg, depth, alsoDual = True):
    """Every quipu-with-relations algebra a bounded walk out of `pathAlg` reaches.

    Returns a dict from `certificate` to the mutation path that got there.  What
    comes back is a set of algebras *proved* to be in the class of `pathAlg`,
    every one of them a member of the family `search` enumerates -- so the
    intersection of the two says how much of the family a class accounts for, and
    the difference says what the polynomial matched without a path to show for
    it.

    The walk is run from the opposite algebra too, for the reason
    `search.linesReachedFrom` gives, and what it reaches is carried back through
    the opposite so the certificates are comparable.
    """
    from . import search

    found = {}

    def visit(quiver, mutationVertices):
        key = certificate(quiver)
        if key is not None and (key not in found or len(mutationVertices) < len(found[key])):
            found[key] = list(mutationVertices)

    def visitDual(quiver, mutationVertices):
        visit(pathAlgebra.dualPathAlgebra(quiver), mutationVertices)

    search.mutationSearchDepthFirst(pathAlg, depth, [], 'quipus', printOutput = False,
                                    visitor = visit)
    if alsoDual:
        search.mutationSearchDepthFirst(
            pathAlgebra.dualPathAlgebra(pathAlg), depth, [], 'quipus',
            printOutput = False, visitor = visitDual)
    return found


def describeCertificate(certificate):
    """A certificate read back as (quipu name, arrows, relations), numbered from 1.

    `reachedQuipuAlgebras` returns certificates rather than algebras, since that
    is what makes them comparable; this turns one back into something to read or
    to draw.
    """
    arrows, relations = certificate
    graph = nx.Graph()
    graph.add_edges_from(arrows)
    parameters = quipuForms.quipuParameters(graph)
    return (
        quipuForms.formatQuipu(parameters),
        tuple("{0}->{1}".format(tail + 1, head + 1) for tail, head in arrows),
        tuple("-".join(str(vertex + 1) for vertex in path) for path in relations),
    )
