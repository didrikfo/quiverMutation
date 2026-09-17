"""Every tree as a hereditary algebra, against every LNA of the same length.

The quipu theorem says an LNA with almost separate relations is derived
equivalent to the path algebra of a quipu, and F-031 asked whether that shape is
forced: it enumerated the trees of maximum degree three that are *not* quipus,
up to order 12, and found none sharing a Coxeter polynomial with any LNA.  This
module asks the same question of **every** tree, degree four and above included,
which is the part F-031 left out.

The mechanics are cheap because a tree carries no relations.  Its path algebra is
hereditary, all orientations of it are derived equivalent by BGP reflection, and
the Cartan matrix is then the reachability matrix of any one orientation.  So
one polynomial per tree, and `networkx.nonisomorphic_trees` enumerates the trees
themselves -- 47 at order 9, 106 at 10, 235 at 11, 551 at 12.

What a match would mean.  The tree's algebra is hereditary, so a match with an
LNA that is in **no quipu class** would be a second hereditary family behind the
Nakayama classification -- the thing the quipu theorem is for a different shape
of tree.  A match with an LNA the theorem already places is a cospectral
coincidence and nothing more, and the report keeps the two apart.
"""

import networkx as nx

from . import coxeterTables
from . import invariants
from . import pathAlgebra
from . import quipuForms


def treesOfOrder(order):
    """Every tree on `order` vertices up to isomorphism, as undirected graphs."""
    if order <= 0:
        return []
    if order <= 2:
        return [nx.path_graph(order)]
    return list(nx.nonisomorphic_trees(order))


def orientedTreeQuiver(graph):
    """A path algebra on the tree, with every edge oriented away from a root.

    Which orientation does not matter -- the orientations of a tree quiver are
    related by BGP reflections, hence all derived equivalent, and they all have
    the same Coxeter polynomial -- so the cheapest one is taken: the BFS tree out
    of the lowest-numbered vertex.
    """
    algebra = pathAlgebra.PathAlgebra()
    algebra.add_vertices_from(sorted(graph.nodes))
    root = min(graph.nodes)
    for source, target in nx.bfs_edges(graph, root):
        algebra.add_arrow(source, target)
    return algebra


def treeCoxeterKey(graph):
    """The Coxeter polynomial of the tree's path algebra, as a coefficient tuple.

    Computed from the reachability matrix of the oriented tree rather than by
    counting paths, which is the same matrix -- in a tree there is at most one
    path between two vertices -- and avoids building the algebra at all.
    """
    nodes = sorted(graph.nodes)
    position = {node: index for index, node in enumerate(nodes)}
    size = len(nodes)
    cartan = [[0] * size for _ in range(size)]
    root = nodes[0]
    parent = {root: None}
    order = [root]
    for source, target in nx.bfs_edges(graph, root):
        parent[target] = source
        order.append(target)
    for node in order:
        cartan[position[node]][position[node]] = 1
        ancestor = parent[node]
        while ancestor is not None:
            cartan[position[node]][position[ancestor]] = 1
            ancestor = parent[ancestor]
    return invariants.coxeterCoefficients(cartan)


def describeTree(graph):
    """A name for the tree: its quipu notation where it has one, else its AHU form."""
    parameters = quipuForms.quipuParameters(graph)
    if parameters is not None:
        return quipuForms.formatQuipu(parameters)
    return quipuForms.canonicalTreeForm(graph)


def treeReport(order):
    """Every tree of the order, matched against the LNAs of that length.

    Returns a list of dicts, one per tree, in the enumeration's order:

    | key | meaning |
    |---|---|
    | `name` | the quipu notation, or the AHU encoding of the tree |
    | `maxDegree` | the tree's maximum degree |
    | `isQuipu` | whether the quipu theorem's shape covers it |
    | `key` | the Coxeter polynomial, as a coefficient tuple |
    | `lnas` | the LNAs of the length sharing that polynomial |
    | `byStatus` | those LNAs grouped by `coxeterTables.lnaStatus` |
    | `cospectralQuipus` | the quipus of the order carrying the same polynomial |

    `cospectralQuipus` is what makes a match readable.  Two path algebras of
    trees are derived equivalent exactly when the trees are isomorphic, so a
    **non-quipu** tree sharing its polynomial with a quipu is cospectral with it
    and derived equivalent to nothing the quipu is derived equivalent to.  Any
    LNA under that polynomial that the moves place in the quipu's class is
    therefore ruled out for this tree on the spot, with no search: the polynomial
    matches and the algebras still cannot be equivalent.
    """
    index = coxeterTables.lnaKeyIndex(order)
    status = coxeterTables.lnaStatus(order)
    quipuIndex = coxeterTables.quipuKeyIndex(order)
    rows = []
    for graph in treesOfOrder(order):
        key = treeCoxeterKey(graph)
        lnas = index.get(key, ())
        byStatus = {}
        for relLengths in lnas:
            byStatus.setdefault(status[relLengths], []).append(relLengths)
        rows.append({
            'name': describeTree(graph),
            'maxDegree': max((degree for _, degree in graph.degree()), default = 0),
            'isQuipu': quipuForms.isQuipuByDegrees(graph),
            'key': key,
            'lnas': lnas,
            'byStatus': {name: tuple(members) for name, members in byStatus.items()},
            'cospectralQuipus': quipuIndex.get(key, ()),
        })
    return rows


def leads(order):
    """The rows worth following: a non-quipu tree against an LNA no quipu holds.

    A match is only a lead when the LNA on the other side is *not* already in a
    quipu class -- status NOT_QUIPU or UNPLACED.  A match against an LNA the
    moves place in a quipu class is not a lead at all but a refutation, by the
    argument in `treeReport`: the tree is cospectral with that quipu and not
    isomorphic to it, so the two hereditary algebras are not derived equivalent
    and neither is the LNA and this tree.

    Empty is the expected answer, and is the result -- F-031 for maximum degree
    three, and this module for the rest.
    """
    found = []
    for row in treeReport(order):
        if row['isQuipu']:
            continue
        open_ = {name: members for name, members in row['byStatus'].items()
                 if name != coxeterTables.QUIPU}
        if open_:
            found.append(dict(row, byStatus = open_))
    return found


def cospectralWithAQuipu(order):
    """Non-quipu trees carrying a polynomial some quipu of the order also carries.

    Not leads -- see `leads` -- but the reason the plain polynomial comparison
    has hits at all, and worth counting: they are the cospectral pairs of F-010
    with one side outside the quipu shape.
    """
    return [row for row in treeReport(order)
            if not row['isQuipu'] and row['cospectralQuipus']]
