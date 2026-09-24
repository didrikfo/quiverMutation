"""Which quiver shapes the mutation walks pass through, counted over every walk.

Three families of non-line quivers have mattered to the classification --
quipus with relations (F-034), commutative squares with a side of two (F-027),
quivers with one parallel pair (H-016) -- and each was found by looking at one
walk by hand.  This asks the question the other way round: walk out of a planned
set of LNAs, record every quiver reached under the relabelling-invariant keys
of `shapeKeys`, and count which shapes are reached from many classes (hubs),
which tie together orbits the classification already joins (bridges), and
which are reached from two classes at all (candidate merges, replayed before
they are believed).

The census half is `walkStart`, a pure function of one start so that
`batch.py atlas` can farm it out; everything about exact keys is left to the
analysis half, which resolves a whole ledger through one `ShapeIndex`.  Spec:
docs/superpowers/specs/2026-09-24-shape-atlas-design.md; hypothesis H-022.
"""

import collections
import copy
import functools
import random

from . import fingerprint
from . import freeMoves
from . import nakayama
from . import pathAlgebra
from . import quipuForms
from . import search
from . import shapeKeys


# -- the census -------------------------------------------------------------

def rowString(row):
    """A relation-length row as the digit string the repo names LNAs by."""
    row = tuple(row)
    if any(arrows > 9 for arrows in row):
        raise ValueError("a relation of ten or more arrows has no digit: {0}".format(row))
    return ''.join(str(arrows) for arrows in row)


def ledgerPath(length, depth, sample, seed):
    """Where a census lives.  Every parameter changes the work, so every one is
    in the name; `--jobs` and the budget do not, so they are not."""
    return "logs/atlas-n{0}-d{1}-s{2}-r{3}.jsonl".format(length, depth, sample, seed)


def walkStart(length, row, depth):
    """Every quiver a depth-`depth` walk out of one LNA reaches, and every step.

    The start is walked, and so is its opposite algebra, with what that reaches
    carried back through the opposite -- exactly as `search.quiversReachedFrom`
    does -- so every node is a quiver the start itself reaches and a step of the
    dual walk is recorded as a *left* mutation, a negative vertex.

    Nodes are deduplicated label-exactly (`shapeKeys.labelId`), keeping the
    shortest path met.  An edge joins the node at a path's prefix to the node
    at the path: the walk is depth-first and visits a node before its children,
    so the prefix is always already recorded.
    """
    start = nakayama.LinearNakayamaAlgebra(length, row)
    nodes = []
    position = {}
    edges = set()

    def record(quiver, path):
        ident = shapeKeys.labelId(quiver)
        if ident not in position:
            position[ident] = len(nodes)
            nodes.append({
                'id': ident,
                'path': list(path),
                'depth': len(path),
                'quiver': shapeKeys.serialise(quiver),
                'features': shapeKeys.features(quiver),
                'buckets': [shapeKeys.bucketOf(quiver, level) for level in (0, 1, 2)],
            })
        elif len(path) < nodes[position[ident]]['depth']:
            nodes[position[ident]]['path'] = list(path)
            nodes[position[ident]]['depth'] = len(path)
        return position[ident]

    for dualised in (False, True):
        atPath = {}

        def visit(quiver, mutationVertices, dualised = dualised, atPath = atPath):
            if dualised:
                quiver = pathAlgebra.dualPathAlgebra(quiver)
            path = [-vertex for vertex in mutationVertices] if dualised else list(mutationVertices)
            here = record(quiver, path)
            atPath[tuple(mutationVertices)] = here
            if mutationVertices:
                parent = atPath.get(tuple(mutationVertices[:-1]))
                if parent is not None:
                    edges.add((parent, path[-1], here))

        begin = pathAlgebra.dualPathAlgebra(start) if dualised else copy.deepcopy(start)
        search.mutationSearchDepthFirst(begin, depth, [], 'atlas', printOutput = False,
                                        visitor = visit, visited = fingerprint.Visited())
    return {
        'start': row,
        'depth': depth,
        'nodes': nodes,
        'edges': [list(edge) for edge in sorted(edges)],
    }


@functools.lru_cache(maxsize = None)
def startTags(length):
    """Every LNA of a length, as its row string, to `(orbit, cls)`.

    The orbit is the move orbit of `freeMoves.derivedOrbits` -- free move, edge
    moves and the double mutation, no rule table, which adds nothing on top of
    those (F-032) -- closed under the relation dual.  The class is the orbit's
    quipu where it holds a row the quipu theorem names, and the orbit itself
    otherwise.  So `cls` is as fine as what is known: two leftover orbits that
    are in truth one class (H-013) are two labels here, and a shape that joins
    them is exactly what the candidate merges are for.
    """
    cover = freeMoves.coverage(length, rules = (), free = True, edges = True, doubles = True)
    parent = {}

    def find(node):
        parent.setdefault(node, node)
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    orbitOf = {}
    for root, members in cover['orbits'].items():
        for member in members:
            orbitOf[member] = root
    for lna in cover['lnas']:
        first, second = find(orbitOf[lna]), find(orbitOf[freeMoves.mirrorRow(length, lna)])
        if first != second:
            parent[first] = second
    groups = collections.defaultdict(list)
    for lna in cover['lnas']:
        groups[find(orbitOf[lna])].append(lna)

    tags = {}
    for members in groups.values():
        orbit = 'orbit:' + rowString(min(members))
        seeded = sorted(member for member in members if member in cover['seeded'])
        if seeded:
            k, m = quipuForms.quipuForAlmostSeparateLNA(length, list(seeded[0]))
            cls = 'quipu:' + quipuForms.formatQuipu(quipuForms.canonicalQuipuParameters(k, m))
        else:
            cls = orbit
        for member in members:
            tags[rowString(member)] = (orbit, cls)
    return tags


def startsFor(length, sample = 0, seed = 0):
    """The rows a census starts from.

    `sample = 0` is every LNA.  Otherwise every row outside a quipu class, and
    up to `sample` rows of each quipu class drawn with a fixed seed -- the
    leftovers are the point, the sample is the contrast.
    """
    tags = startTags(length)
    rows = sorted(tags)
    if not sample:
        return rows
    byClass = collections.defaultdict(list)
    for row in rows:
        byClass[tags[row][1]].append(row)
    chooser = random.Random(seed)
    chosen = []
    for cls in sorted(byClass):
        members = byClass[cls]
        if cls.startswith('orbit:') or len(members) <= sample:
            chosen.extend(members)
        else:
            chosen.extend(sorted(chooser.sample(members, sample)))
    return chosen
