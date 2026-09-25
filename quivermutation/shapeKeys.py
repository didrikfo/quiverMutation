"""Keys for a quiver with relations that do not see its labels.

`search.quiverKey` and `fingerprint.canonicalKey` are exact *with labels*, and
inside one walk that is right: vertex labels do not move under mutation
(F-049).  Across walks it is wrong.  Two LNAs whose walks reach the same quiver
under different vertex labels are never seen to meet, and a count of "which
shapes do the walks pass through" cannot be made at all.  This module is the
relabelling-invariant half, at four levels of detail:

| level | keeps |
|---|---|
| L0 | the underlying undirected multigraph |
| L1 | the quiver, parallel arrows included |
| L2 | the quiver and a relation *skeleton*: each relation's kind and path lengths, as a hyperedge on the vertices it runs through |
| L3 | the algebra as presented, up to relabelling, parallel-arrow naming and the sign gauge |

Each level is a labelled undirected graph -- one node per vertex, one per arrow,
one per relation -- whose edge labels carry the direction (`tail`, `head`), so
an isomorphism of the graph is an isomorphism of the quiver.  The
Weisfeiler-Lehman hash of that graph is a **bucket**, not a key: WL cannot tell
some non-isomorphic graphs apart (two triangles and a hexagon are the classic
pair).  `ShapeIndex` resolves a bucket into exact keys by testing isomorphism
against the representatives already in it.

L3 is decided by `fingerprint.canonicalKey` after relabelling through each
isomorphism of the L2 graphs, so it inherits that key's one conservatism: two
presentations of one ideal by genuinely different generators are two keys.  A
missed match costs a meeting; there are no false ones.  Spec:
docs/superpowers/specs/2026-09-24-shape-atlas-design.md.
"""

import functools
import hashlib
import json
from fractions import Fraction

import networkx as nx
from networkx.algorithms import isomorphism

from . import arrowPaths
from . import fingerprint
from . import procedure
from . import quipuForms


LEVELS = (0, 1, 2, 3)

#: Rounds of WL refinement.  The graphs are small (tens of nodes) and four
#: rounds reach across any relation the walks produce at these lengths.
WL_ITERATIONS = 4


def relationKind(relation):
    """'zero' for a monomial, 'comm' for two terms, 'other' for more."""
    terms = len(relation)
    if terms == 1:
        return 'zero'
    if terms == 2:
        return 'comm'
    return 'other'


def structureGraph(pathAlg, level):
    """The labelled graph a level compares, as an undirected `nx.Graph`.

    Level 3 has no graph of its own: it is compared through the level-2 graph's
    isomorphisms (see `ShapeIndex`), so asking for it gives the level-2 graph.
    """
    level = min(level, 2)
    graph = nx.Graph()
    for vertex in pathAlg.quiver.nodes:
        graph.add_node(('v', vertex), label = 'v')
    for tail, head, key in pathAlg.quiver.edges(keys = True):
        arrow = ('a', tail, head, key)
        graph.add_node(arrow, label = 'a')
        # At L0 both ends look the same, which is what forgets the orientation.
        graph.add_edge(('v', tail), arrow, label = 'e' if level == 0 else 'tail')
        graph.add_edge(arrow, ('v', head), label = 'e' if level == 0 else 'head')
    if level < 2:
        return graph
    for index, relation in enumerate(procedure.relationsFrom(pathAlg)):
        if not relation:
            continue
        paths = sorted(relation)
        lengths = sorted(len(path) for path in paths)
        node = ('r', index)
        graph.add_node(node, label = '{0}:{1}'.format(
            relationKind(relation), ','.join(map(str, lengths))))
        source = arrowPaths.pathSource(paths[0])
        target = arrowPaths.pathTarget(paths[0])
        graph.add_edge(node, ('v', source), label = 'start')
        graph.add_edge(node, ('v', target), label = 'end')
        interior = set()
        for path in paths:
            interior.update(arrowPaths.pathVertices(path)[1:-1])
        for vertex in interior:
            graph.add_edge(node, ('v', vertex), label = 'through')
    return graph


def bucketOf(pathAlg, level):
    """The WL hash of a level's graph: equal for isomorphic shapes, and a bucket
    rather than a key because it can also be equal for non-isomorphic ones."""
    return nx.weisfeiler_lehman_graph_hash(
        structureGraph(pathAlg, level), node_attr = 'label', edge_attr = 'label',
        iterations = WL_ITERATIONS)


def _coefficientText(value):
    return str(Fraction(value))


def _coefficientValue(text):
    value = Fraction(text)
    return int(value) if value.denominator == 1 else value


def serialise(pathAlg):
    """The algebra as JSON-safe data: vertices, arrows with keys, arrow relations."""
    return {
        'vertices': sorted(pathAlg.quiver.nodes),
        'arrows': [list(arrow) for arrow in arrowPaths.arrowsOf(pathAlg.quiver)],
        'relations': [
            [[[list(arrow) for arrow in path], _coefficientText(coefficient)]
             for path, coefficient in sorted(relation.items())]
            for relation in procedure.relationsFrom(pathAlg)
        ],
    }


def deserialise(data):
    """`serialise` read back.  Integral coefficients come back as `int`, since
    `fingerprint.digest` hashes `repr` and `Fraction(1)` does not repr as `1`."""
    quiver = nx.MultiDiGraph()
    quiver.add_nodes_from(data['vertices'])
    for tail, head, key in data['arrows']:
        quiver.add_edge(tail, head, key = key)
    relations = [
        {tuple(tuple(arrow) for arrow in path): _coefficientValue(coefficient)
         for path, coefficient in relation}
        for relation in data['relations']
    ]
    return procedure.toPathAlgebra(quiver, relations)


def relabel(pathAlg, mapping):
    """The same algebra with vertex `v` renamed `mapping[v]`, arrow keys kept."""
    quiver = nx.MultiDiGraph()
    quiver.add_nodes_from(mapping[vertex] for vertex in pathAlg.quiver.nodes)
    for tail, head, key in pathAlg.quiver.edges(keys = True):
        quiver.add_edge(mapping[tail], mapping[head], key = key)
    relations = [
        {tuple((mapping[tail], mapping[head], key) for tail, head, key in path): coefficient
         for path, coefficient in relation.items()}
        for relation in procedure.relationsFrom(pathAlg)
    ]
    return procedure.toPathAlgebra(quiver, relations)


def labelId(pathAlg):
    """A label-exact identity: the digest of `fingerprint.canonicalKey`.

    Where the key is refused (a bundle structure past its cap, never yet seen),
    the serialised presentation is hashed instead -- the conservative direction,
    since two presentations of one algebra then get two ids.
    """
    key = fingerprint.canonicalKey(pathAlg)
    if key is None:
        raw = json.dumps(serialise(pathAlg), sort_keys = True)
        return 'raw:' + hashlib.blake2b(raw.encode('utf-8'), digest_size = 16).hexdigest()
    return '{0:032x}'.format(fingerprint.digest(key))


def squareSides(pathAlg):
    """The commutative squares, as 'short x long' in arrows, or None.

    A square is a commutativity relation whose two paths share only their ends
    -- F-027's shape, which the walks pass through on the way back to a line.
    Several are joined by '+', sorted.
    """
    found = []
    for relation in procedure.relationsFrom(pathAlg):
        if len(relation) != 2:
            continue
        first, second = sorted(relation)
        inner = set(arrowPaths.pathVertices(first)[1:-1])
        if inner & set(arrowPaths.pathVertices(second)[1:-1]):
            continue
        short, long = sorted((len(first), len(second)))
        found.append('{0}x{1}'.format(short, long))
    return '+'.join(sorted(found)) if found else None


def features(pathAlg):
    """Cheap numbers about one quiver, computed where it is recorded."""
    quiver = pathAlg.quiver
    vertexCount = quiver.number_of_nodes()
    arrowCount = quiver.number_of_edges()
    undirected = quiver.to_undirected()
    components = nx.number_connected_components(undirected) if vertexCount else 0
    simple = nx.Graph(undirected)
    bundles = fingerprint.arrowBundles(quiver)
    isTree = (not bundles and vertexCount > 0 and nx.is_tree(simple))
    isLine = isTree and all(quiver.in_degree(v) <= 1 and quiver.out_degree(v) <= 1
                            for v in quiver.nodes)
    isQuipu = isTree and quipuForms.isQuipuByDegrees(simple)
    relationCount = len(pathAlg.rels)
    cords = sum(1 for _vertex, degree in simple.degree() if degree == 3)
    return {
        'vertices': vertexCount,
        'arrows': arrowCount,
        'relations': relationCount,
        'bundles': len(bundles),
        'cycleRank': arrowCount - vertexCount + components,
        'isLine': isLine,
        'isTree': isTree,
        'isQuipu': isQuipu,
        'quipu': quipuForms.formatQuipu(quipuForms.quipuParameters(simple)) if isQuipu else None,
        # H-017's defect, relations less cords, where the shape is a quipu.
        'defect': relationCount - cords if isQuipu else None,
        'square': squareSides(pathAlg),
    }


#: The most isomorphisms of two L2 graphs `ShapeIndex` will try when asking
#: whether one of them carries one algebra onto the other.  A shape with a
#: larger automorphism group than this is not expected at these lengths; if it
#: happens the comparison says "different", which costs a meeting and never
#: invents one, and `ShapeIndex.capHits` counts it so it is not silent.
ISOMORPHISM_CAP = 5000

_NODE_MATCH = isomorphism.categorical_node_match('label', None)
_EDGE_MATCH = isomorphism.categorical_edge_match('label', None)


def _matcher(first, second):
    return isomorphism.GraphMatcher(first, second, node_match = _NODE_MATCH,
                                    edge_match = _EDGE_MATCH)


#: How many representatives' structure graphs `ShapeIndex` keeps built.  Most
#: buckets hold one shape and are never compared again, so a small cache of the
#: recently compared ones serves nearly every comparison.
GRAPH_CACHE = 4096


def compactForm(pathAlg):
    """The algebra as one compact JSON string: what `ShapeIndex` keeps of a
    representative.  A built `nx.Graph` is several KB; at n = 9 there are about
    130 000 representatives per level, which is how the n = 9 analysis ran out
    of memory at 6.1 GB (task 8a)."""
    return json.dumps(serialise(pathAlg), separators = (',', ':'))


@functools.lru_cache(maxsize = GRAPH_CACHE)
def _representativeGraph(text, level):
    # `serialise` keeps vertex names, arrow keys and every relation's support,
    # so the graph built back from it is the representative's own graph, node
    # for node -- which L3 needs, since it relabels through the node names.
    return structureGraph(deserialise(json.loads(text)), level)


class ShapeIndex:
    """Exact keys at every level, one representative per shape per bucket.

    `keyOf` hashes the quiver into its bucket and compares it with the bucket's
    representatives: by graph isomorphism at L0 to L2, and at L3 by relabelling
    the algebra through each isomorphism of the L2 graphs and comparing
    `fingerprint.canonicalKey`.  An isomorphism of algebras preserves every
    relation's kind, lengths and support, so it is always among the L2
    isomorphisms, and trying those is complete up to the cap.

    A representative is kept as its `compactForm` string (and at L3 its
    canonical key), and its graph is built again only when a later quiver
    lands in the same bucket.  Vertex names must therefore survive JSON, as
    the integers every census uses do.

    Keys are only comparable within one index: they number representatives in
    the order they were met.
    """

    def __init__(self, isomorphismCap = ISOMORPHISM_CAP):
        self.isomorphismCap = isomorphismCap
        self.capHits = 0
        self._representatives = {level: {} for level in LEVELS}

    def keyOf(self, pathAlg, level, bucket = None):
        bucket = bucketOf(pathAlg, level) if bucket is None else bucket
        representatives = self._representatives[level].setdefault(bucket, [])
        # The quiver's own graph is only needed when there is something to
        # compare it with; a first-met bucket, the common case, never builds it.
        graph = structureGraph(pathAlg, level) if representatives else None
        if level < 3:
            for position, representative in enumerate(representatives):
                if _matcher(graph, _representativeGraph(representative, level)).is_isomorphic():
                    return self._key(level, bucket, position)
            representatives.append(compactForm(pathAlg))
            return self._key(level, bucket, len(representatives) - 1)
        canonical = fingerprint.canonicalKey(pathAlg)
        for position, (representative, representativeKey) in enumerate(representatives):
            if self._sameAlgebra(pathAlg, canonical, graph, representative, representativeKey):
                return self._key(level, bucket, position)
        representatives.append((compactForm(pathAlg), canonical))
        return self._key(level, bucket, len(representatives) - 1)

    def _sameAlgebra(self, pathAlg, canonical, graph, representative, representativeKey):
        if canonical is None or representativeKey is None:
            return False
        if canonical == representativeKey:
            return True
        representativeGraph = _representativeGraph(representative, 3)
        for count, mapping in enumerate(_matcher(graph, representativeGraph).isomorphisms_iter()):
            if count >= self.isomorphismCap:
                self.capHits += 1
                return False
            renaming = {node[1]: image[1] for node, image in mapping.items() if node[0] == 'v'}
            if fingerprint.canonicalKey(relabel(pathAlg, renaming)) == representativeKey:
                return True
        return False

    @staticmethod
    def _key(level, bucket, position):
        return 'L{0}:{1}:{2}'.format(level, bucket, position)


def describe(pathAlg):
    """One line a person can read: the arrows, then each relation."""
    arrows = ' '.join('{0}->{1}'.format(tail, head)
                      for tail, head, _key in arrowPaths.arrowsOf(pathAlg.quiver))
    relations = []
    for relation in procedure.relationsFrom(pathAlg):
        paths = ['-'.join(map(str, arrowPaths.pathVertices(path))) for path in sorted(relation)]
        relations.append('{0} {1}'.format(relationKind(relation), ' = '.join(paths)))
    return '{0} | {1}'.format(arrows, '; '.join(relations) or 'no relations')
