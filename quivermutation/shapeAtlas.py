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
import json
import random

import polars as pl

from . import fingerprint
from . import freeMoves
from . import invariants
from . import mutation
from . import nakayama
from . import pathAlgebra
from . import quipuForms
from . import quipuRelations
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


# -- the analysis -------------------------------------------------------------

def iterLedger(path):
    """A ledger's records one line at a time, skipping a line that will not
    parse as `jobs.Ledger.records` does.  The n = 9 ledger is 332 MB and parsed
    whole it no longer fits beside the analysis (task 8a)."""
    try:
        handle = open(path)
    except FileNotFoundError:
        return
    with handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                yield json.loads(line)
            except ValueError:
                continue


def resolve(records, length, index = None):
    """A ledger as three tables, every distinct quiver keyed at every level.

    One `ShapeIndex` over the whole ledger, so keys are comparable across
    starts -- which is the point -- and each label-exact quiver is resolved once
    however many starts reached it.

    `records` may be a list or a stream (`iterLedger`): each record is read
    once, and of each distinct quiver only what the tables need is kept -- its
    JSON as a string, its buckets and its features -- so memory grows with the
    distinct quivers and the visits, never with the ledger's text.
    """
    tags = startTags(length)
    index = shapeKeys.ShapeIndex() if index is None else index
    # One column list per table column rather than a dict per row; an id met
    # again is replaced by the string first stored for it, so every visit and
    # edge of it points at one object.
    known = {}
    nodeColumns = {'id': [], 'quiver': []}
    buckets = []
    featureColumns = None
    visitColumns = {name: [] for name in ('start', 'orbit', 'cls', 'id', 'depth', 'path')}
    edgeColumns = {name: [] for name in ('start', 'parent', 'vertex', 'child')}
    for record in records:
        start = record['unit']
        result = record['result']
        orbit, cls = tags[start]
        ids = []
        for node in result['nodes']:
            seen = len(known)
            ident = known.setdefault(node['id'], node['id'])
            ids.append(ident)
            if len(known) > seen:
                nodeColumns['id'].append(ident)
                nodeColumns['quiver'].append(json.dumps(node['quiver']))
                buckets.append(tuple(node['buckets']))
                if featureColumns is None:
                    featureColumns = {name: [] for name in node['features']}
                for name, column in featureColumns.items():
                    column.append(node['features'][name])
            visitColumns['start'].append(start)
            visitColumns['orbit'].append(orbit)
            visitColumns['cls'].append(cls)
            visitColumns['id'].append(ident)
            visitColumns['depth'].append(node['depth'])
            visitColumns['path'].append(','.join(str(step) for step in node['path']))
        for parent, vertex, child in result['edges']:
            edgeColumns['start'].append(start)
            edgeColumns['parent'].append(ids[parent])
            edgeColumns['vertex'].append(vertex)
            edgeColumns['child'].append(ids[child])
    del known

    # The visits and edges are finished: make them frames, and let the lists go
    # before the index starts to grow.
    visitSchema = {'start': pl.Utf8, 'orbit': pl.Utf8, 'cls': pl.Utf8, 'id': pl.Utf8,
                   'depth': pl.Int64, 'path': pl.Utf8}
    edgeSchema = {'start': pl.Utf8, 'parent': pl.Utf8, 'vertex': pl.Int64, 'child': pl.Utf8}
    visits = pl.DataFrame(visitColumns, schema = visitSchema)
    del visitColumns
    edges = pl.DataFrame(edgeColumns, schema = edgeSchema)
    del edgeColumns

    keyColumns = {'key{0}'.format(level): [] for level in shapeKeys.LEVELS}
    for text, nodeBuckets in zip(nodeColumns['quiver'], buckets):
        algebra = shapeKeys.deserialise(json.loads(text))
        for level in shapeKeys.LEVELS:
            keyColumns['key{0}'.format(level)].append(
                index.keyOf(algebra, level, bucket = nodeBuckets[min(level, 2)]))
    del buckets
    columns = [pl.Series(name, values, dtype = pl.Utf8)
               for name, values in list(nodeColumns.items()) + list(keyColumns.items())]
    # A feature column is typed from all its values, as the row-wise frame it
    # replaces was (`infer_schema_length = None`): some are None for most rows.
    columns += [pl.Series(name, values, strict = False)
                for name, values in (featureColumns or {}).items()]
    return {
        'nodes': pl.DataFrame(columns),
        'visits': visits,
        'edges': edges,
    }


def writeTables(tables, stem):
    for name, table in tables.items():
        table.write_parquet('{0}.{1}.parquet'.format(stem, name))


def returnRates(tables, level):
    """Per shape, the share of the starts reaching it whose walk goes on from it
    to a line **other than the start**, without passing back through the start.

    The start is excluded on purpose: the dual walk records the left mutations
    back to it, so every node would otherwise "return" trivially.
    """
    key = 'key{0}'.format(level)
    nodes = tables['nodes']
    keyOf = dict(zip(nodes['id'].to_list(), nodes[key].to_list()))
    lines = set(nodes.filter(pl.col('isLine'))['id'].to_list())
    roots = dict(tables['visits'].filter(pl.col('depth') == 0).select('start', 'id').iter_rows())
    parentsByStart = collections.defaultdict(lambda: collections.defaultdict(set))
    for start, parent, child in tables['edges'].select('start', 'parent', 'child').iter_rows():
        parentsByStart[start][child].add(parent)
    seen = collections.Counter()
    returning = collections.Counter()
    for start, ids in tables['visits'].group_by('start').agg(pl.col('id')).iter_rows():
        root = roots[start]
        parents = parentsByStart[start]
        good = set()
        frontier = [ident for ident in set(ids) if ident in lines and ident != root]
        while frontier:
            node = frontier.pop()
            for parent in parents[node]:
                if parent != root and parent not in good:
                    good.add(parent)
                    frontier.append(parent)
        for shape in {keyOf[ident] for ident in ids}:
            seen[shape] += 1
        for shape in {keyOf[ident] for ident in good}:
            returning[shape] += 1
    return pl.DataFrame({key: list(seen), 'returnRate': [returning[s] / seen[s] for s in seen]},
                        schema = {key: pl.Utf8, 'returnRate': pl.Float64})


def shapeMeasures(tables, level):
    """Per shape at one level: how widely, how early, how mixed, and whether it
    leads back to a line.  Sorted by classes reached, then starts."""
    key = 'key{0}'.format(level)
    visits = tables['visits'].join(tables['nodes'].select('id', key), on = 'id')
    perStart = visits.group_by(key, 'start', 'orbit', 'cls').agg(pl.col('depth').min())
    base = perStart.group_by(key).agg(
        pl.col('start').n_unique().alias('starts'),
        pl.col('orbit').n_unique().alias('orbits'),
        pl.col('cls').n_unique().alias('classes'),
        pl.col('cls').str.starts_with('orbit:').mean().alias('leftoverShare'),
        pl.col('depth').min().alias('firstDepth'),
        pl.col('depth').median().alias('medianDepth'),
    )
    shares = perStart.group_by(key, 'cls').len().with_columns(
        (pl.col('len') / pl.col('len').sum().over(key)).alias('p'))
    mixing = shares.group_by(key).agg(
        (-(pl.col('p') * pl.col('p').log(2))).sum().abs().alias('mixing'))
    # A bridge: some class reaches the shape from two or more of its orbits, so
    # the shape sits where the classification already joined them.
    bridges = perStart.group_by(key, 'cls').agg(
        pl.col('orbit').n_unique().alias('orbitsInClass')).group_by(key).agg(
        (pl.col('orbitsInClass').max() > 1).alias('bridge'))
    isLine = tables['nodes'].group_by(key).agg(pl.col('isLine').any())
    return (base.join(mixing, on = key).join(bridges, on = key)
            .join(returnRates(tables, level), on = key, how = 'left')
            .join(isLine, on = key)
            .with_columns(pl.col('returnRate').fill_null(0.0))
            .sort(['classes', 'starts'], descending = True))


def transitions(tables, level = 2, maxLength = 4, top = 30, keep = 100):
    """The commonest cycles line -> S1 -> ... -> line through the shape graph.

    Every line is one token, `LINE`, whatever its relations; every other node is
    its shape at `level`.  A step is counted once per distinct label-exact pair
    of quivers.  A cycle's weight is its weakest step, and only the `keep` most
    connected shapes are searched through, which is what bounds the search.
    These are rule *templates* in F-027's sense: open, walk, close.
    """
    key = 'key{0}'.format(level)
    nodes = tables['nodes']
    token = {ident: ('LINE' if line else shape)
             for ident, shape, line in nodes.select('id', key, 'isLine').iter_rows()}
    counts = collections.Counter(
        (token[parent], token[child])
        for parent, child in tables['edges'].select('parent', 'child').unique().iter_rows())
    out = collections.defaultdict(dict)
    weight = collections.Counter()
    for (first, second), count in counts.items():
        if first == second:
            continue
        out[first][second] = count
        weight[first] += count
        weight[second] += count
    kept = {shape for shape, _count in weight.most_common(keep)} | {'LINE'}
    cycles = []

    def extend(path, bottleneck):
        for following, count in out[path[-1]].items():
            if following not in kept:
                continue
            narrowest = min(bottleneck, count)
            if following == 'LINE':
                if len(path) >= 2:
                    cycles.append((narrowest, path[1:]))
            elif following not in path and len(path) < maxLength:
                extend(path + [following], narrowest)

    extend(['LINE'], float('inf'))
    cycles.sort(key = lambda cycle: (-cycle[0], len(cycle[1]), cycle[1]))
    return cycles[:top]


def describeKey(tables, level, key):
    """One representative of a shape, readably."""
    quiver = tables['nodes'].filter(pl.col('key{0}'.format(level)) == key)['quiver'][0]
    return shapeKeys.describe(shapeKeys.deserialise(json.loads(quiver)))


def _short(key):
    level, bucket, position = key.split(':')
    return '{0}:{1}:{2}'.format(level, bucket[:8], position)


def report(tables, level, top, out):
    """What a ledger shows, as text: counts, hubs, bridges, leftover hubs, cycles."""
    nodes = tables['nodes']
    visits = tables['visits']
    print("{0} starts, {1} classes, {2} distinct quivers with their labels".format(
        visits['start'].n_unique(), visits['cls'].n_unique(), nodes.height), file = out)
    print("shapes: " + ", ".join("L{0} {1}".format(l, nodes['key{0}'.format(l)].n_unique())
                                 for l in shapeKeys.LEVELS), file = out)
    measures = shapeMeasures(tables, level)
    key = 'key{0}'.format(level)
    columns = ['starts', 'orbits', 'classes', 'firstDepth', 'returnRate', 'leftoverShare']

    def table(title, frame):
        print("\n{0} (L{1})".format(title, level), file = out)
        for row in frame.head(top).iter_rows(named = True):
            print("  {0:<24} {1}".format(_short(row[key]), "  ".join(
                "{0}={1:.2f}".format(c, row[c]) if isinstance(row[c], float)
                else "{0}={1}".format(c, row[c]) for c in columns)), file = out)
            print("      " + describeKey(tables, level, row[key]), file = out)

    nonLines = measures.filter(~pl.col('isLine'))
    table("hubs", nonLines)
    table("bridges", nonLines.filter(pl.col('bridge')))
    table("hubs among leftovers", nonLines.filter(pl.col('leftoverShare') > 0)
          .sort('leftoverShare', 'starts', descending = True))
    print("\ncycles through a line (L{0})".format(level), file = out)
    for count, shapes in transitions(tables, level, top = top):
        print("  {0:>6}  line -> {1} -> line".format(
            count, " -> ".join(_short(shape) for shape in shapes)), file = out)


# -- candidate merges -----------------------------------------------------------

def _side(row):
    return {'start': row['start'], 'orbit': row['orbit'], 'cls': row['cls'], 'id': row['id'],
            'path': [int(step) for step in row['path'].split(',') if step]}


def candidateMerges(tables):
    """Pairs of classes whose walks reach one L3 shape, shortest paths first.

    Each is only a candidate until `replay` has re-run it: the key could be
    wrong, and a key bug would show here first.  One candidate per pair of
    classes, the cheapest meeting kept.
    """
    visits = tables['visits'].join(tables['nodes'].select('id', 'key3'), on = 'id')
    shared = (visits.group_by('key3').agg(pl.col('cls').n_unique().alias('classes'))
              .filter(pl.col('classes') > 1).select('key3'))
    rows = visits.join(shared, on = 'key3').sort('key3', 'depth', 'start')
    firstPerClass = collections.defaultdict(dict)
    for row in rows.iter_rows(named = True):
        firstPerClass[row['key3']].setdefault(row['cls'], row)
    best = {}
    for key, byClass in firstPerClass.items():
        classes = sorted(byClass)
        for one, other in zip(classes, classes[1:]):
            cost = byClass[one]['depth'] + byClass[other]['depth']
            if (one, other) not in best or cost < best[(one, other)]['cost']:
                best[(one, other)] = {'key': key, 'cost': cost,
                                      'first': _side(byClass[one]),
                                      'second': _side(byClass[other])}
    return sorted(best.values(),
                  key = lambda c: (c['cost'], c['first']['cls'], c['second']['cls']))


def replay(length, candidate):
    """Re-run both paths step by step and check they end at one algebra.

    Every step must be admissible (for a left step, on the opposite algebra, as
    the dual walk took it) and keep the start's Coxeter polynomial, as the
    search's guard requires (F-038) -- a cyclic end with no polynomial is let
    through, as the search lets it through.  Then the two ends must have the
    same L3 key in a fresh index, and the two starts the same polynomial.
    """
    ends = []
    for side in (candidate['first'], candidate['second']):
        start = nakayama.LinearNakayamaAlgebra(length, side['start'])
        baseKey = invariants.coxeterKey(start)
        algebra = copy.deepcopy(start)
        for step in side['path']:
            checked = algebra if step > 0 else pathAlgebra.dualPathAlgebra(algebra)
            if not mutation.mutationIsPossibleAtVertex(checked, abs(step)):
                return {'ok': False, 'reason': 'step {0} of {1} not admissible'.format(
                    step, side['start'])}
            algebra = mutation.quiverMutationAtVertices(algebra, [step])
            moved = search._coxeterKeyOrNone(algebra)
            if moved is not None and moved != baseKey:
                return {'ok': False, 'reason': 'step {0} of {1} moved the Coxeter polynomial'
                        .format(step, side['start'])}
        ends.append((algebra, baseKey))
    (first, firstKey), (second, secondKey) = ends
    if firstKey != secondKey:
        return {'ok': False, 'reason': 'the two starts have different Coxeter polynomials'}
    index = shapeKeys.ShapeIndex()
    if index.keyOf(first, 3) != index.keyOf(second, 3):
        return {'ok': False, 'reason': 'the replayed ends are not isomorphic'}
    return {'ok': True, 'reason': 'replayed'}


# -- validation -----------------------------------------------------------------

#: H-014's hub for the n = 9 leftovers: the line on eight vertices with one
#: pendant vertex at the second, carrying relations.  Named the way
#: `shapeKeys.features` names a quipu.
N9_HUB = 'P^(6)_(1,1)'


def _squaresVerdict(shortSides, topSquare):
    """Pure verdict for the squares check, split out of `validate` so the rule
    can be tested without running a census.  F-027 holds when short side 2 is
    both the top square's own short side and strictly the most common one
    among all returning squares; H-022's stronger "only ever 2" clause was
    refuted by E-053, so a minority of other short sides (e.g. 3xk) is fine."""
    if not shortSides:
        return False, None
    majorityShortSide = shortSides.most_common(1)[0][0]
    isStrictMajority = shortSides[2] > 0 and all(
        shortSides[2] > count for side, count in shortSides.items() if side != 2)
    ok = (majorityShortSide == 2 and isStrictMajority
          and topSquare is not None and topSquare.startswith('2x'))
    return ok, majorityShortSide


def validate(tables, length, records, coverageSample = 40):
    """The four checks of H-022.  A failure of the first two says the instrument
    is wrong; of the last two, that the keys are.

    `records` is only read by check 4, for a sample of starts: a list of
    records, or the ledger's path, which is then re-read for just the sample."""
    result = {}
    nodes = tables['nodes']

    # 1. F-027's dominance of short side 2 among the squares of the shapes that
    # lead back to a line.  H-022 had also predicted no short side of 3 at
    # all, but the n = 8 census (E-053) found 36 genuine 3xk returning
    # squares (e.g. 405000 by [3, 2]) alongside 957 with a short side of 2, so
    # that clause is refuted; F-027's "never a 3" was only ever measured on
    # steps inside rule-table rules, and over-generalised from there.
    measures = shapeMeasures(tables, 2)
    squareOf = nodes.group_by('key2').agg(pl.col('square').first())
    returning = (measures.filter(~pl.col('isLine') & (pl.col('returnRate') > 0))
                 .join(squareOf, on = 'key2').sort('starts', descending = True))
    shortSides = collections.Counter()
    for square in returning['square'].drop_nulls().to_list():
        for part in square.split('+'):
            shortSides[int(part.split('x')[0])] += 1
    withSquare = returning.filter(pl.col('square').is_not_null())
    topSquare = withSquare['square'][0] if withSquare.height else None
    ok, majorityShortSide = _squaresVerdict(shortSides, topSquare)
    result['squares'] = {
        'returningShapes': returning.height,
        'shortSides': dict(sorted(shortSides.items())),
        'topSquare': topSquare,
        'majorityShortSide': majorityShortSide,
        'ok': ok,
    }

    # 2. H-014's quipu hub for the leftovers, at n = 9 only.
    if length == 9:
        hubIds = nodes.filter((pl.col('quipu') == N9_HUB) & (pl.col('relations') > 0))['id']
        leftoverStarts = (tables['visits'].filter(pl.col('cls').str.starts_with('orbit:')
                                                  & pl.col('id').is_in(hubIds.to_list()))
                          ['start'].n_unique())
        result['quipuHub'] = {'name': N9_HUB, 'leftoverStarts': leftoverStarts,
                              'ok': leftoverStarts >= 7}
    else:
        result['quipuHub'] = 'skipped'

    # 3. On quipus with monomial relations the L3 key and the certificate of
    #    `quipuRelations` are two exact routes to one answer.
    keyToCertificates = collections.defaultdict(set)
    certificateToKeys = collections.defaultdict(set)
    for key, quiver in nodes.filter(pl.col('isQuipu')).select('key3', 'quiver').iter_rows():
        certificate = quipuRelations.certificate(shapeKeys.deserialise(json.loads(quiver)))
        if certificate is None:
            continue
        keyToCertificates[key].add(certificate)
        certificateToKeys[certificate].add(key)
    violations = ([key for key, found in keyToCertificates.items() if len(found) > 1]
                  + [str(c) for c, found in certificateToKeys.items() if len(found) > 1])
    result['certificates'] = {'checked': len(keyToCertificates), 'violations': violations,
                              'ok': not violations}

    # 4. The census reaches everything the label-exact search does, so every
    #    meeting that search can find is a shared id here, hence a shared key.
    sampled = _coverageSample(records, coverageSample)
    missing = 0
    for record in sampled:
        got = {search.quiverKey(shapeKeys.deserialise(node['quiver']))
               for node in record['result']['nodes']}
        start = nakayama.LinearNakayamaAlgebra(length, record['unit'])
        depth = min(2, record['result']['depth'])
        expected = set(search.quiversReachedFrom(start, depth, alsoDual = True)) - {None}
        missing += len(expected - got)
    result['coverage'] = {'startsChecked': len(sampled), 'missing': missing,
                          'ok': missing == 0}
    return result


def _coverageSample(records, coverageSample):
    """Every `step`-th record in order of start, as check 4 samples them.

    `records` is a list, or a ledger path: then the ledger is read twice, once
    for the start names and once for just the sampled records, so the whole
    ledger is never held at once.
    """
    if not isinstance(records, str):
        ordered = sorted(records, key = lambda record: record['unit'])
        return ordered[::max(1, len(ordered) // coverageSample)]
    units = sorted(record['unit'] for record in iterLedger(records))
    wanted = collections.Counter(units[::max(1, len(units) // coverageSample)])
    found = []
    for record in iterLedger(records):
        if wanted[record['unit']] > 0:
            wanted[record['unit']] -= 1
            found.append(record)
    return sorted(found, key = lambda record: record['unit'])
