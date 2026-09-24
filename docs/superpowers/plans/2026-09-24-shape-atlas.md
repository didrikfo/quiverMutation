# Shape Atlas Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Record every quiver the depth-4 mutation walks reach out of the LNAs of lengths 8 to 10, key each one at four levels up to relabelling, and read off which shapes are hubs and bridges between classes — with any shape shared by two classes replayed as a candidate merge.

**Architecture:** `quivermutation/shapeKeys.py` turns a quiver with relations into labelled structure graphs, Weisfeiler–Lehman buckets and exact keys (L0–L3), with a `ShapeIndex` that resolves buckets into keys by exact isomorphism. `quivermutation/shapeAtlas.py` holds the census walk (one start per unit), the start tags, and the analysis over a ledger (tables, measures, transitions, candidate merges, replay, validation). `batch.py atlas` is the census front door, `atlas.py` the analysis front door, and `quivermutation/atlasPage.py` draws the top shapes as a static page.

**Tech Stack:** Python ≥ 3.10, networkx 3.6 (WL hash, VF2 `GraphMatcher`), polars 1.44 (tables, parquet), pytest. No new dependencies.

**Spec:** `docs/superpowers/specs/2026-09-24-shape-atlas-design.md`

## Global Constraints

- Python and tests run in the project's Linux venv under WSL. From Git Bash, every command below is prefixed with
  `MSYS_NO_PATHCONV=1 wsl.exe -d Ubuntu --cd /mnt/c/Users/didri/kode/quiverMutation -- ` and uses `.venv/bin/python`. This prefix is written `$WSL` below, e.g. `$WSL .venv/bin/python -m pytest tests/test_shape_keys.py -q`. (Define it once per shell: `WSL="MSYS_NO_PATHCONV=1 wsl.exe -d Ubuntu --cd /mnt/c/Users/didri/kode/quiverMutation --"` and run as `eval $WSL .venv/bin/python ...`, or type the prefix out.)
- Dependencies: only what `pyproject.toml` already lists (matplotlib, networkx, numpy, polars, sympy). **No scipy** — so no `kamada_kawai_layout`; use `spring_layout`.
- Repo style: camelCase function and variable names, module docstrings that say *why*, comments explaining reasoning, research identifiers (F-/H-/E-/R-) cited where a design choice comes from one.
- New tests must be fast (seconds). Anything over ~10 s gets `@pytest.mark.slow`.
- Every ledger lives under `logs/` (gitignored); no data files are committed.
- Research harness conventions (`research/README.md`): new entries go at the **top** of their file, dated, with an identifier; nothing is deleted. Next free identifiers: **H-022**, **E-053**, **F-054**.
- Commit messages follow the repo's style (a plain sentence saying what changed and why) and end with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Work on branch `claude/shape-atlas` (already checked out; the spec is committed there).

## File map

| file | responsibility |
|---|---|
| `quivermutation/shapeKeys.py` (new) | structure graphs per level, WL buckets, exact keys (`ShapeIndex`), serialise/deserialise, relabel, label-exact id, cheap features, a readable description |
| `quivermutation/shapeAtlas.py` (new) | census walk of one start, start tags and start selection, ledger path; analysis: tables, measures, return rates, transitions, candidate merges, replay, validation, report |
| `quivermutation/atlasPage.py` (new) | a static HTML page drawing shapes as SVG |
| `batch.py` (modify) | `AtlasTask`, registered in `TASKS`; `atlas.py` listed in `ELSEWHERE` |
| `atlas.py` (new) | command line over a ledger: resolve, write parquet, report, replay, validate, page |
| `pyproject.toml` (modify) | add `atlas` to `py-modules` |
| `tests/test_shape_keys.py` (new) | keys |
| `tests/test_shape_atlas.py` (new) | census, task, analysis, replay, validation, page |
| `research/HYPOTHESES.md`, `GLOSSARY.md`, `README.md` (modify) | H-022, terms, a README section |
| `research/EXPERIMENTS.md` (modify, Task 8) | E-053, the runs |

## Where the plan departs from the spec, and why

- **No Coxeter polynomial per node.** The Coxeter guard stops the walk from taking any step that changes the polynomial (F-038), so every node of one start's walk has the start's polynomial. Storing it per node would add nothing. `replay` still compares the two starts' polynomials.
- **The page does not reuse `classpage`'s drawing.** That code is browser JavaScript that draws lines and quipus from their names. Atlas shapes have no names, so `atlasPage` draws them in Python from the quiver.
- **How the return rate is defined.** It is "goes on to a line other than the start, without passing back through the start". The dual walk records the left mutations back to the start, so without that exclusion every node would count as returning.
- **L3 buckets are the L2 buckets.** An isomorphism of algebras is always an isomorphism of the L2 graphs, so the L2 hash is a valid bucket for L3.

---

### Task 1: Structure graphs, buckets, serialisation, relabelling, features

**Files:**
- Create: `quivermutation/shapeKeys.py`
- Test: `tests/test_shape_keys.py`

**Interfaces:**
- Consumes (existing): `procedure.relationsFrom(pathAlg) -> list[dict[path, coefficient]]` where a path is a tuple of `(tail, head, key)` arrows and coefficients are `int` or `Fraction`; `procedure.toPathAlgebra(quiver, relations) -> PathAlgebra`; `arrowPaths.arrowsOf(quiver)`, `arrowPaths.pathVertices(path)`; `fingerprint.canonicalKey(pathAlg)`, `fingerprint.digest(key) -> int`, `fingerprint.arrowBundles(quiver)`; `quipuForms.isQuipuByDegrees(graph)`, `quipuForms.quipuParameters(graph)`, `quipuForms.formatQuipu(params)`.
- Produces:
  - `LEVELS = (0, 1, 2, 3)`, `WL_ITERATIONS = 4`
  - `relationKind(relation: dict) -> str` — `'zero' | 'comm' | 'other'`
  - `structureGraph(pathAlg, level: int) -> nx.Graph` — level clamped to 0..2
  - `bucketOf(pathAlg, level: int) -> str` — WL hash; level 3 uses the level-2 graph
  - `serialise(pathAlg) -> dict` (JSON-safe), `deserialise(data: dict) -> PathAlgebra`
  - `relabel(pathAlg, mapping: dict) -> PathAlgebra`
  - `labelId(pathAlg) -> str` — label-exact identity (32 hex chars, or `'raw:'`+32 hex)
  - `squareSides(pathAlg) -> str | None` — e.g. `'2x4'`, several joined by `'+'`
  - `features(pathAlg) -> dict` with keys `vertices, arrows, relations, bundles, cycleRank, isLine, isTree, isQuipu, quipu, defect, square`
  - `describe(pathAlg) -> str`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_shape_keys.py`:

```python
"""Keys for a quiver with relations that do not see its labels.

The atlas (spec docs/superpowers/specs/2026-09-24-shape-atlas-design.md) counts
shapes across walks out of *different* LNAs, where the same quiver turns up
under different vertex labels.  Everything here is about the one property that
makes that count mean anything: a key changes when the shape changes and not
when the labels do.
"""

import json
import random
from fractions import Fraction

import networkx as nx
import pytest

from quivermutation import fingerprint
from quivermutation import mutation
from quivermutation import nakayama as nk
from quivermutation import procedure
from quivermutation import shapeKeys as sk


def _square():
    """F-027's example: A_10 with one relation on 3..8, right mutation at 3.

    The result has a zero relation 2-4-3 and a commutativity square between
    4-3-8 (two arrows) and 4-5-6-7-8 (four).
    """
    return mutation.quiverMutationAtVertices(nk.LinearNakayamaAlgebra(10, "00500000"), [3])


def _quiver(arrows):
    quiver = nx.MultiDiGraph()
    for tail, head, key in arrows:
        quiver.add_edge(tail, head, key = key)
    return quiver


def _permuted(pathAlg, seed):
    vertices = sorted(pathAlg.quiver.nodes)
    images = vertices[:]
    random.Random(seed).shuffle(images)
    return sk.relabel(pathAlg, dict(zip(vertices, images)))


@pytest.mark.parametrize("seed", range(5))
def test_relabelling_leaves_every_bucket_unchanged(seed):
    algebra = _square()
    moved = _permuted(algebra, seed)
    for level in sk.LEVELS:
        assert sk.bucketOf(algebra, level) == sk.bucketOf(moved, level)


def test_relabel_is_a_relabelling_and_not_a_different_algebra():
    algebra = _square()
    identity = {vertex: vertex for vertex in algebra.quiver.nodes}
    assert fingerprint.canonicalKey(sk.relabel(algebra, identity)) == fingerprint.canonicalKey(algebra)
    moved = _permuted(algebra, 1)
    assert sorted(moved.quiver.nodes) == sorted(algebra.quiver.nodes)
    assert moved.quiver.number_of_edges() == algebra.quiver.number_of_edges()


def test_zero_and_commutativity_agree_at_L1_and_differ_at_L2():
    arrows = [(1, 2, 0), (2, 4, 0), (1, 3, 0), (3, 4, 0)]
    comm = procedure.toPathAlgebra(_quiver(arrows), [
        {((1, 2, 0), (2, 4, 0)): 1, ((1, 3, 0), (3, 4, 0)): -1}])
    zero = procedure.toPathAlgebra(_quiver(arrows), [
        {((1, 2, 0), (2, 4, 0)): 1}, {((1, 3, 0), (3, 4, 0)): 1}])
    assert sk.bucketOf(comm, 1) == sk.bucketOf(zero, 1)
    assert sk.bucketOf(comm, 2) != sk.bucketOf(zero, 2)


def test_serialise_round_trips_through_json():
    algebra = _square()
    back = sk.deserialise(json.loads(json.dumps(sk.serialise(algebra))))
    assert fingerprint.canonicalKey(back) == fingerprint.canonicalKey(algebra)
    assert sk.labelId(back) == sk.labelId(algebra)


def test_deserialise_keeps_integral_coefficients_integral():
    """`fingerprint.digest` hashes `repr`, and `repr(Fraction(1)) != repr(1)`."""
    back = sk.deserialise(json.loads(json.dumps(sk.serialise(_square()))))
    for relation in procedure.relationsFrom(back):
        for coefficient in relation.values():
            assert not (isinstance(coefficient, Fraction) and coefficient.denominator == 1)


def test_features_of_a_line():
    found = sk.features(nk.LinearNakayamaAlgebra(6, "3000"))
    assert found['isLine'] and found['isTree'] and found['isQuipu']
    assert found['cycleRank'] == 0
    assert found['relations'] == 1
    assert found['square'] is None


def test_features_of_the_square():
    found = sk.features(_square())
    assert not found['isLine'] and not found['isTree'] and not found['isQuipu']
    assert found['cycleRank'] == 1
    assert found['square'] == '2x4'
    assert found['defect'] is None


def test_features_of_a_quipu_count_its_defect():
    tree = procedure.toPathAlgebra(_quiver([(1, 2, 0), (2, 3, 0), (2, 4, 0)]), [])
    found = sk.features(tree)
    assert found['isQuipu'] and not found['isLine']
    assert found['defect'] == -1          # no relations, one cord
    assert found['quipu'].startswith('P^(')


def test_describe_names_arrows_and_relations():
    text = sk.describe(_square())
    assert '4->3' in text and 'comm' in text and 'zero' in text
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `$WSL .venv/bin/python -m pytest tests/test_shape_keys.py -q`
Expected: FAIL at collection with `ImportError: cannot import name 'shapeKeys'`.

- [ ] **Step 3: Write the implementation**

Create `quivermutation/shapeKeys.py`:

```python
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


def describe(pathAlg):
    """One line a person can read: the arrows, then each relation."""
    arrows = ' '.join('{0}->{1}'.format(tail, head)
                      for tail, head, _key in arrowPaths.arrowsOf(pathAlg.quiver))
    relations = []
    for relation in procedure.relationsFrom(pathAlg):
        paths = ['-'.join(map(str, arrowPaths.pathVertices(path))) for path in sorted(relation)]
        relations.append('{0} {1}'.format(relationKind(relation), ' = '.join(paths)))
    return '{0} | {1}'.format(arrows, '; '.join(relations) or 'no relations')
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `$WSL .venv/bin/python -m pytest tests/test_shape_keys.py -q`
Expected: all pass. If `test_features_of_the_square` reports a different `square` string, print `sk.describe(_square())` and check it against F-027's example (`(4,3,8) = (4,5,6,7,8)`) before changing the test.

- [ ] **Step 5: Commit**

```bash
git add quivermutation/shapeKeys.py tests/test_shape_keys.py
git commit -m "Key a quiver with relations by its shape, at four levels, without its labels

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Exact keys — `ShapeIndex`

**Files:**
- Modify: `quivermutation/shapeKeys.py` (append)
- Test: `tests/test_shape_keys.py` (append)

**Interfaces:**
- Consumes: everything Task 1 produced.
- Produces:
  - `ISOMORPHISM_CAP = 5000`
  - `class ShapeIndex` with `__init__(self, isomorphismCap = ISOMORPHISM_CAP)`, `keyOf(self, pathAlg, level, bucket = None) -> str` returning `'L{level}:{bucket}:{index}'`, and attribute `capHits: int` (how many L3 comparisons stopped at the cap — each is a possibly missed match).

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_shape_keys.py`:

```python
@pytest.mark.parametrize("seed", range(5))
def test_relabelling_leaves_every_key_unchanged(seed):
    algebra = _square()
    moved = _permuted(algebra, seed)
    index = sk.ShapeIndex()
    for level in sk.LEVELS:
        assert index.keyOf(algebra, level) == index.keyOf(moved, level)


def test_a_WL_collision_still_gets_two_keys():
    """Two triangles and a hexagon: WL cannot separate them, VF2 must."""
    triangles = procedure.toPathAlgebra(_quiver([
        (1, 2, 0), (2, 3, 0), (1, 3, 0), (4, 5, 0), (5, 6, 0), (4, 6, 0)]), [])
    hexagon = procedure.toPathAlgebra(_quiver([
        (1, 2, 0), (2, 3, 0), (3, 4, 0), (4, 5, 0), (5, 6, 0), (1, 6, 0)]), [])
    assert sk.bucketOf(triangles, 0) == sk.bucketOf(hexagon, 0)
    index = sk.ShapeIndex()
    assert index.keyOf(triangles, 0) != index.keyOf(hexagon, 0)


def test_zero_and_commutativity_differ_as_keys_at_L2_and_L3():
    arrows = [(1, 2, 0), (2, 4, 0), (1, 3, 0), (3, 4, 0)]
    comm = procedure.toPathAlgebra(_quiver(arrows), [
        {((1, 2, 0), (2, 4, 0)): 1, ((1, 3, 0), (3, 4, 0)): -1}])
    zero = procedure.toPathAlgebra(_quiver(arrows), [
        {((1, 2, 0), (2, 4, 0)): 1}, {((1, 3, 0), (3, 4, 0)): 1}])
    index = sk.ShapeIndex()
    assert index.keyOf(comm, 1) == index.keyOf(zero, 1)
    assert index.keyOf(comm, 2) != index.keyOf(zero, 2)
    assert index.keyOf(comm, 3) != index.keyOf(zero, 3)


def _kronecker(relations):
    return procedure.toPathAlgebra(_quiver([(1, 2, 0), (1, 2, 1), (2, 3, 0)]), relations)


def test_parallel_naming_and_sign_gauge_agree_at_L3():
    first = _kronecker([{((1, 2, 0), (2, 3, 0)): 1}])
    second = _kronecker([{((1, 2, 1), (2, 3, 0)): 1}])
    plus = _kronecker([{((1, 2, 0), (2, 3, 0)): 1, ((1, 2, 1), (2, 3, 0)): 1}])
    minus = _kronecker([{((1, 2, 0), (2, 3, 0)): 1, ((1, 2, 1), (2, 3, 0)): -1}])
    index = sk.ShapeIndex()
    assert index.keyOf(first, 3) == index.keyOf(second, 3)
    assert index.keyOf(plus, 3) == index.keyOf(minus, 3)
    assert index.keyOf(first, 3) != index.keyOf(plus, 3)


def test_two_different_lnas_are_two_keys_at_L3_and_one_at_L0():
    first = nk.LinearNakayamaAlgebra(6, "3000")
    second = nk.LinearNakayamaAlgebra(6, "0300")
    index = sk.ShapeIndex()
    assert index.keyOf(first, 0) == index.keyOf(second, 0)
    assert index.keyOf(first, 3) != index.keyOf(second, 3)


def test_a_key_carries_its_level():
    assert sk.ShapeIndex().keyOf(_square(), 2).startswith('L2:')
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `$WSL .venv/bin/python -m pytest tests/test_shape_keys.py -q`
Expected: the new tests FAIL with `AttributeError: module 'quivermutation.shapeKeys' has no attribute 'ShapeIndex'`.

- [ ] **Step 3: Write the implementation**

Append to `quivermutation/shapeKeys.py`:

```python
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


class ShapeIndex:
    """Exact keys at every level, one representative per shape per bucket.

    `keyOf` hashes the quiver into its bucket and compares it with the bucket's
    representatives: by graph isomorphism at L0 to L2, and at L3 by relabelling
    the algebra through each isomorphism of the L2 graphs and comparing
    `fingerprint.canonicalKey`.  An isomorphism of algebras preserves every
    relation's kind, lengths and support, so it is always among the L2
    isomorphisms, and trying those is complete up to the cap.

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
        graph = structureGraph(pathAlg, level)
        if level < 3:
            for position, representative in enumerate(representatives):
                if _matcher(graph, representative).is_isomorphic():
                    return self._key(level, bucket, position)
            representatives.append(graph)
            return self._key(level, bucket, len(representatives) - 1)
        canonical = fingerprint.canonicalKey(pathAlg)
        for position, (representativeGraph, representativeKey) in enumerate(representatives):
            if self._sameAlgebra(pathAlg, canonical, graph,
                                 representativeGraph, representativeKey):
                return self._key(level, bucket, position)
        representatives.append((graph, canonical))
        return self._key(level, bucket, len(representatives) - 1)

    def _sameAlgebra(self, pathAlg, canonical, graph, representativeGraph, representativeKey):
        if canonical is None or representativeKey is None:
            return False
        if canonical == representativeKey:
            return True
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
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `$WSL .venv/bin/python -m pytest tests/test_shape_keys.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add quivermutation/shapeKeys.py tests/test_shape_keys.py
git commit -m "Resolve shape buckets into exact keys, by isomorphism and the canonical key

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: The census — walk one start, tag the starts, `batch.py atlas`

**Files:**
- Create: `quivermutation/shapeAtlas.py`
- Modify: `batch.py` (add `AtlasTask`; register it in `TASKS` at `batch.py:1018`)
- Test: `tests/test_shape_atlas.py`

**Interfaces:**
- Consumes: `shapeKeys.labelId`, `serialise`, `features`, `bucketOf`; existing `search.mutationSearchDepthFirst(pathAlg, depth, [], name, printOutput = False, visitor = fn, visited = fingerprint.Visited())` where `fn(pathAlg, mutationVertices)` is called at every node before its children; `pathAlgebra.dualPathAlgebra`; `freeMoves.coverage(length, rules = (), free = True, edges = True, doubles = True)` returning a dict with `'lnas'`, `'orbits'` (root → list of row tuples), `'seeded'` (set of row tuples); `freeMoves.mirrorRow(length, row)`; `quipuForms.quipuForAlmostSeparateLNA(length, list(row)) -> (k, m)`, `quipuForms.canonicalQuipuParameters(k, m)`, `quipuForms.formatQuipu`.
- Produces:
  - `rowString(row) -> str`
  - `ledgerPath(length, depth, sample, seed) -> str` = `'logs/atlas-n{length}-d{depth}-s{sample}-r{seed}.jsonl'`
  - `walkStart(length, row: str, depth: int) -> dict` with keys `start, depth, nodes, edges`; each node is `{'id', 'path', 'depth', 'quiver', 'features', 'buckets'}` (`buckets` = `[b0, b1, b2]`); each edge is `[parentIndex, signedVertex, childIndex]` (negative vertex = left mutation, from the dual walk). Node 0 is the start itself.
  - `startTags(length) -> dict[str, tuple[str, str]]` row string → `(orbit, cls)`; `orbit` is `'orbit:<least member>'` of the dual-closed move orbit, `cls` is `'quipu:P^(..)_(..)'` when the orbit holds a seeded row and otherwise equals `orbit`.
  - `startsFor(length, sample = 0, seed = 0) -> list[str]`
  - `batch.AtlasTask` (`name = 'atlas'`)

- [ ] **Step 1: Write the failing tests**

Create `tests/test_shape_atlas.py`:

```python
"""The shape atlas: the census walk, the task, and the analysis over a ledger.

Checked at n = 5 and 6, where every LNA is in a quipu class and a depth-2 walk
takes a fraction of a second, so the instrument can be pinned against things
already known before it is run where nothing is.  Spec:
docs/superpowers/specs/2026-09-24-shape-atlas-design.md.
"""

import argparse
import io
import json

import pytest

import batch
from quivermutation import jobs
from quivermutation import nakayama as nk
from quivermutation import search
from quivermutation import shapeAtlas as sa
from quivermutation import shapeKeys as sk


def test_a_walk_records_everything_the_label_exact_search_reaches():
    record = sa.walkStart(5, "300", 2)
    got = {search.quiverKey(sk.deserialise(node['quiver'])) for node in record['nodes']}
    expected = set(search.quiversReachedFrom(nk.LinearNakayamaAlgebra(5, "300"), 2,
                                             alsoDual = True))
    assert expected - {None} <= got


def test_the_first_node_is_the_start_and_every_edge_joins_recorded_nodes():
    record = sa.walkStart(5, "300", 2)
    first = record['nodes'][0]
    assert first['path'] == [] and first['depth'] == 0 and first['features']['isLine']
    size = len(record['nodes'])
    assert record['edges']
    for parent, vertex, child in record['edges']:
        assert 0 <= parent < size and 0 <= child < size and vertex != 0
    assert len({node['id'] for node in record['nodes']}) == size


def test_a_walk_is_json():
    record = sa.walkStart(5, "300", 2)
    assert json.loads(json.dumps(record)) == record


def test_every_n6_lna_is_tagged_with_a_quipu_class_and_mirrors_share_an_orbit():
    tags = sa.startTags(6)
    assert len(tags) == 42                                   # Catalan(5)
    assert all(cls.startswith('quipu:P^(') for _orbit, cls in tags.values())
    for row in tags:
        mirror = sa.rowString(batch.fm.mirrorRow(6, tuple(int(c) for c in row)))
        assert tags[row][0] == tags[mirror][0]


def test_sampling_keeps_every_leftover_and_one_per_quipu_class():
    tags = sa.startTags(9)
    chosen = sa.startsFor(9, sample = 1, seed = 0)
    leftovers = sorted(row for row, (_orbit, cls) in tags.items() if cls.startswith('orbit:'))
    assert len(leftovers) == 9                               # F-032, n = 9
    assert set(leftovers) <= set(chosen)
    quipuClasses = {cls for _orbit, cls in tags.values() if cls.startswith('quipu:')}
    assert len(chosen) == len(leftovers) + len(quipuClasses)
    assert sa.startsFor(9, sample = 1, seed = 0) == chosen   # deterministic


def test_the_task_writes_a_ledger_and_resumes_from_it(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    task = batch.TASKS['atlas']
    args = argparse.Namespace(length = 5, depth = 1, sample = 0, seed = 0)
    assert jobs.runTask(task, args, out = io.StringIO()) == 0
    records = jobs.Ledger(task.ledgerPath(args)).records()
    assert len(records) == 14                                # Catalan(4)
    assert records[0]['result']['nodes'][0]['depth'] == 0
    again = io.StringIO()
    assert jobs.runTask(task, args, out = again) == 0
    assert "nothing to do" in again.getvalue()
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `$WSL .venv/bin/python -m pytest tests/test_shape_atlas.py -q`
Expected: FAIL at collection with `ImportError: cannot import name 'shapeAtlas'`.

- [ ] **Step 3: Write `quivermutation/shapeAtlas.py` (census half)**

```python
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
```

- [ ] **Step 4: Add `AtlasTask` to `batch.py`**

Add `from quivermutation import shapeAtlas` to the imports at the top of `batch.py` (after `from quivermutation import sampling`). Insert this class directly above the line `TASKS = {task.name: task for task in [SampleTask(), CoresTask()]}` (currently `batch.py:1018`), and change that line to include it:

```python
class AtlasTask(jobs.Task):
    """Walk out of LNAs and record every quiver shape the walks pass through.

    One unit is one starting LNA: a depth-`--depth` walk out of it and out of
    its opposite, every quiver reached recorded label-exactly with its WL
    buckets, features and shortest path, and every step between them.  The
    ledger is read by `python atlas.py`, which resolves the buckets into exact
    keys up to relabelling and counts hubs, bridges and candidate merges.

    `--sample 0` walks every LNA of the length, which is what n = 8 and 9 want.
    `--sample N` walks every LNA outside a quipu class and N of each quipu class,
    which is what n = 10 wants.  Spec:
    docs/superpowers/specs/2026-09-24-shape-atlas-design.md; H-022.
    """

    name = 'atlas'
    help = "record every quiver shape the walks out of a length's LNAs reach"

    def addArguments(self, parser):
        parser.add_argument("length", type = int, help = "the line length")
        parser.add_argument("--depth", type = int, default = 4,
                            help = "mutations to walk out of each start (default 4)")
        parser.add_argument("--sample", type = int, default = 0,
                            help = "rows per quipu class; 0 walks every LNA (default 0)")
        parser.add_argument("--seed", type = int, default = 0,
                            help = "the sample's seed (default 0)")

    def ledgerPath(self, args):
        return shapeAtlas.ledgerPath(args.length, args.depth, args.sample, args.seed)

    def units(self, args):
        return shapeAtlas.startsFor(args.length, args.sample, args.seed)

    def run(self, unit, args):
        return shapeAtlas.walkStart(args.length, unit, args.depth)

    def summarise(self, records, args, out = sys.stdout):
        results = [record['result'] for record in records]
        if not results:
            print("nothing in the ledger yet", file = out)
            return
        recorded = sum(len(result['nodes']) for result in results)
        distinct = len({node['id'] for result in results for node in result['nodes']})
        print("n = {0}, depth {1}: {2} starts walked".format(
            args.length, args.depth, len(results)), file = out)
        print("  {0} quivers recorded, {1} distinct with their labels".format(
            recorded, distinct), file = out)
        print("  read it with: python atlas.py {0} --depth {1} --sample {2} --seed {3}".format(
            args.length, args.depth, args.sample, args.seed), file = out)


TASKS = {task.name: task for task in [SampleTask(), CoresTask(), AtlasTask()]}
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `$WSL .venv/bin/python -m pytest tests/test_shape_atlas.py -q`
Expected: all pass. `test_sampling_keeps_every_leftover_and_one_per_quipu_class` runs the n = 9 coverage and may take several seconds; if it takes over 10 s, mark it `@pytest.mark.slow`.

- [ ] **Step 6: Commit**

```bash
git add quivermutation/shapeAtlas.py batch.py tests/test_shape_atlas.py
git commit -m "Walk out of each LNA and record every quiver reached, as a batch task

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: The analysis — tables, measures, return rates, transitions, report

**Files:**
- Modify: `quivermutation/shapeAtlas.py` (append)
- Test: `tests/test_shape_atlas.py` (append)

**Interfaces:**
- Consumes: Task 3's ledger records (`{'unit': row, 'result': walkStart(...)}`), `startTags`, `shapeKeys.ShapeIndex`, `shapeKeys.deserialise`, `shapeKeys.describe`.
- Produces:
  - `resolve(records, length, index = None) -> dict[str, pl.DataFrame]` with tables:
    - `nodes`: `id, key0, key1, key2, key3, quiver` (JSON string) + every feature column
    - `visits`: `start, orbit, cls, id, depth, path` (path as `'1,-3,2'`)
    - `edges`: `start, parent, vertex, child` (ids)
  - `writeTables(tables, stem)` → `{stem}.{name}.parquet`
  - `returnRates(tables, level) -> pl.DataFrame[key{level}, returnRate]`
  - `shapeMeasures(tables, level) -> pl.DataFrame` with `key{level}, starts, orbits, classes, leftoverShare, firstDepth, medianDepth, mixing, bridge, returnRate, isLine`, sorted by `classes, starts` descending
  - `transitions(tables, level = 2, maxLength = 4, top = 30, keep = 100) -> list[tuple[int, list[str]]]` — `(bottleneck count, [shape keys])` for cycles `LINE → S1 → … → LINE`
  - `describeKey(tables, level, key) -> str`
  - `report(tables, level, top, out)`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_shape_atlas.py`:

```python
import polars as pl


@pytest.fixture(scope = "module")
def atlas5():
    """Every LNA of length 5, walked to depth 2, resolved."""
    records = [{'unit': row, 'result': sa.walkStart(5, row, 2)} for row in sa.startsFor(5)]
    return {'records': records, 'tables': sa.resolve(records, 5)}


def test_resolve_gives_one_row_per_id_and_keys_at_every_level(atlas5):
    nodes = atlas5['tables']['nodes']
    assert nodes['id'].n_unique() == nodes.height
    for level in (0, 1, 2, 3):
        assert nodes['key{0}'.format(level)].null_count() == 0
    visits = atlas5['tables']['visits']
    assert set(visits['id'].to_list()) <= set(nodes['id'].to_list())


def test_coarser_levels_never_split_what_finer_levels_join(atlas5):
    nodes = atlas5['tables']['nodes']
    for finer, coarser in ((3, 2), (2, 1), (1, 0)):
        grouped = nodes.group_by('key{0}'.format(finer)).agg(
            pl.col('key{0}'.format(coarser)).n_unique().alias('n'))
        assert grouped['n'].max() == 1


def test_every_start_has_exactly_one_depth_zero_visit_and_it_is_a_line(atlas5):
    tables = atlas5['tables']
    roots = tables['visits'].filter(pl.col('depth') == 0).join(
        tables['nodes'].select('id', 'isLine'), on = 'id')
    assert roots['start'].n_unique() == roots.height == 14
    assert roots['isLine'].all()


def test_a_line_is_reached_from_its_own_class_only(atlas5):
    measures = sa.shapeMeasures(atlas5['tables'], 3)
    lines = measures.filter(pl.col('isLine'))
    assert lines.height > 0
    assert lines['classes'].max() == 1


def test_measures_are_in_range(atlas5):
    measures = sa.shapeMeasures(atlas5['tables'], 2)
    assert measures['returnRate'].min() >= 0 and measures['returnRate'].max() <= 1
    assert measures['starts'].min() >= 1
    assert (measures['classes'] <= measures['orbits']).all()


def test_transitions_start_and_end_at_a_line(atlas5):
    cycles = sa.transitions(atlas5['tables'], level = 2)
    assert cycles
    for count, shapes in cycles:
        assert count >= 1 and 1 <= len(shapes) <= 3 and 'LINE' not in shapes


def test_report_prints_the_sections(atlas5):
    out = io.StringIO()
    sa.report(atlas5['tables'], 2, 5, out)
    text = out.getvalue()
    assert 'hubs' in text and 'bridges' in text and 'line -> ' in text
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `$WSL .venv/bin/python -m pytest tests/test_shape_atlas.py -q`
Expected: the new tests FAIL with `AttributeError: module 'quivermutation.shapeAtlas' has no attribute 'resolve'`.

- [ ] **Step 3: Write the analysis half**

Append to `quivermutation/shapeAtlas.py` (and add `import json` and `import polars as pl` to its imports):

```python
# -- the analysis -------------------------------------------------------------

def resolve(records, length, index = None):
    """A ledger as three tables, every distinct quiver keyed at every level.

    One `ShapeIndex` over the whole ledger, so keys are comparable across
    starts -- which is the point -- and each label-exact quiver is resolved once
    however many starts reached it.
    """
    tags = startTags(length)
    index = shapeKeys.ShapeIndex() if index is None else index
    nodes = {}
    visits = []
    edges = []
    for record in records:
        start = record['unit']
        result = record['result']
        orbit, cls = tags[start]
        ids = [node['id'] for node in result['nodes']]
        for node in result['nodes']:
            nodes.setdefault(node['id'], node)
            visits.append({'start': start, 'orbit': orbit, 'cls': cls, 'id': node['id'],
                           'depth': node['depth'],
                           'path': ','.join(str(step) for step in node['path'])})
        for parent, vertex, child in result['edges']:
            edges.append({'start': start, 'parent': ids[parent], 'vertex': vertex,
                          'child': ids[child]})
    rows = []
    for ident, node in nodes.items():
        algebra = shapeKeys.deserialise(node['quiver'])
        row = {'id': ident, 'quiver': json.dumps(node['quiver'])}
        for level in shapeKeys.LEVELS:
            row['key{0}'.format(level)] = index.keyOf(algebra, level,
                                                      bucket = node['buckets'][min(level, 2)])
        row.update(node['features'])
        rows.append(row)
    edgeSchema = {'start': pl.Utf8, 'parent': pl.Utf8, 'vertex': pl.Int64, 'child': pl.Utf8}
    return {
        'nodes': pl.DataFrame(rows, infer_schema_length = None),
        'visits': pl.DataFrame(visits, infer_schema_length = None),
        'edges': pl.DataFrame(edges, schema = edgeSchema),
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
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `$WSL .venv/bin/python -m pytest tests/test_shape_atlas.py -q`
Expected: all pass. If `test_a_line_is_reached_from_its_own_class_only` fails, **stop**: either the keys join two different LNAs (a key bug) or the Coxeter-guarded walk crossed classes (a search bug, R-012 territory). Print the offending key's `describeKey` and the starts reaching it before changing anything.

- [ ] **Step 5: Commit**

```bash
git add quivermutation/shapeAtlas.py tests/test_shape_atlas.py
git commit -m "Resolve an atlas ledger into tables, and count hubs, bridges and cycles through a line

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Candidate merges, replay, validation

**Files:**
- Modify: `quivermutation/shapeAtlas.py` (append)
- Test: `tests/test_shape_atlas.py` (append)

**Interfaces:**
- Consumes: Task 4's tables; existing `mutation.mutationIsPossibleAtVertex(pathAlg, vertex)`, `mutation.quiverMutationAtVertices(pathAlg, [step])` (negative step = left mutation, reduces after), `search._coxeterKeyOrNone(pathAlg)`, `invariants.coxeterKey(pathAlg)`, `quipuRelations.certificate(pathAlg)` (exact iso certificate of a quipu with monomial relations, or `None`), `search.quiversReachedFrom(pathAlg, depth, alsoDual = True)`, `search.quiverKey(pathAlg)`.
- Produces:
  - `candidateMerges(tables) -> list[dict]`, each `{'key', 'cost', 'first': side, 'second': side}` with `side = {'start', 'orbit', 'cls', 'id', 'path': list[int]}`
  - `_side(row: dict) -> dict` (builds a side from a visits row)
  - `replay(length, candidate) -> dict` with `ok: bool`, `reason: str`
  - `validate(tables, length, records, coverageSample = 40) -> dict` with sections `squares`, `quipuHub`, `certificates`, `coverage`, each carrying an `ok` flag (or `'skipped'`)

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_shape_atlas.py`:

```python
def test_no_shape_is_shared_by_two_classes_at_n5(atlas5):
    """Every class at n = 5 is known and distinct; a shared L3 key would be a bug."""
    assert sa.candidateMerges(atlas5['tables']) == []


def _sharedNonLineKey(tables):
    visits = tables['visits'].filter(pl.col('depth') > 0).join(
        tables['nodes'].select('id', 'key3', 'isLine'), on = 'id').filter(~pl.col('isLine'))
    shared = visits.group_by('key3').agg(pl.col('start').n_unique().alias('n')).filter(
        pl.col('n') > 1).sort('key3')
    key = shared['key3'][0]
    rows = (visits.filter(pl.col('key3') == key).sort('depth', 'start')
            .unique('start', keep = 'first', maintain_order = True).head(2))
    return key, [sa._side(row) for row in rows.iter_rows(named = True)]


def test_replay_confirms_a_real_shared_key(atlas5):
    key, (first, second) = _sharedNonLineKey(atlas5['tables'])
    result = sa.replay(5, {'key': key, 'first': first, 'second': second})
    assert result['ok'], result


def test_replay_refuses_a_false_one(atlas5):
    key, (first, second) = _sharedNonLineKey(atlas5['tables'])
    result = sa.replay(5, {'key': key, 'first': first, 'second': dict(second, path = [])})
    assert not result['ok']


def test_validation_at_n5(atlas5):
    result = sa.validate(atlas5['tables'], 5, atlas5['records'], coverageSample = 5)
    assert result['certificates']['ok'], result['certificates']
    assert result['coverage']['ok'], result['coverage']
    assert result['quipuHub'] == 'skipped'
    assert 'shortSides' in result['squares']
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `$WSL .venv/bin/python -m pytest tests/test_shape_atlas.py -q`
Expected: the new tests FAIL with `AttributeError: ... has no attribute 'candidateMerges'`.

- [ ] **Step 3: Write the implementation**

Append to `quivermutation/shapeAtlas.py` (add `from . import invariants`, `from . import mutation` and `from . import quipuRelations` to its imports):

```python
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


def validate(tables, length, records, coverageSample = 40):
    """The four checks of H-022.  A failure of the first two says the instrument
    is wrong; of the last two, that the keys are."""
    result = {}
    nodes = tables['nodes']

    # 1. F-027's squares among the shapes that lead back to a line.
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
    result['squares'] = {
        'returningShapes': returning.height,
        'shortSides': dict(sorted(shortSides.items())),
        'topSquare': topSquare,
        'ok': bool(shortSides) and set(shortSides) == {2}
              and topSquare is not None and topSquare.startswith('2x'),
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
    ordered = sorted(records, key = lambda record: record['unit'])
    step = max(1, len(ordered) // coverageSample)
    missing = 0
    for record in ordered[::step]:
        got = {search.quiverKey(shapeKeys.deserialise(node['quiver']))
               for node in record['result']['nodes']}
        start = nakayama.LinearNakayamaAlgebra(length, record['unit'])
        depth = min(2, record['result']['depth'])
        expected = set(search.quiversReachedFrom(start, depth, alsoDual = True)) - {None}
        missing += len(expected - got)
    result['coverage'] = {'startsChecked': len(ordered[::step]), 'missing': missing,
                          'ok': missing == 0}
    return result
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `$WSL .venv/bin/python -m pytest tests/test_shape_atlas.py -q`
Expected: all pass. If `test_replay_confirms_a_real_shared_key` fails with "not isomorphic" on a path with a negative step, the likely cause is that reducing on the opposite algebra (as the dual walk did) and reducing on the algebra itself (as `quiverMutationAtVertices` does) give different presentations. Confirm that by printing `shapeKeys.describe` of both ends. If so, change `replay` to take left steps as `dualPathAlgebra(quiverMutationAtVertices(dualPathAlgebra(algebra), [abs(step)]))`, which is how the walk took them, and re-run.

- [ ] **Step 5: Commit**

```bash
git add quivermutation/shapeAtlas.py tests/test_shape_atlas.py
git commit -m "Find shapes two classes share, replay them, and check the atlas against what is known

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: The page and `atlas.py`

**Files:**
- Create: `quivermutation/atlasPage.py`, `atlas.py`
- Modify: `pyproject.toml` (`py-modules`), `batch.py` (`ELSEWHERE`)
- Test: `tests/test_shape_atlas.py` (append)

**Interfaces:**
- Consumes: `shapeAtlas.resolve, writeTables, report, shapeMeasures, candidateMerges, replay, validate, ledgerPath`; `jobs.Ledger`; `shapeKeys.deserialise, describe`.
- Produces:
  - `atlasPage.drawQuiver(data: dict, size = 240) -> str` (an `<svg>` element)
  - `atlasPage.render(title: str, sections: list[tuple[str, list[dict]]]) -> str` (a full HTML document); each entry dict has `heading`, `quiver` (serialised dict), `lines` (list of str)
  - `atlasPage.sectionsFrom(tables, level, top) -> list[tuple[str, list[dict]]]`
  - `atlas.main(argv = None) -> int`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_shape_atlas.py`:

```python
import atlas
from quivermutation import atlasPage


def test_a_drawing_has_a_node_per_vertex_and_a_path_per_arrow():
    data = sk.serialise(sk.deserialise(json.loads(json.dumps(
        sk.serialise(nk.LinearNakayamaAlgebra(5, "300"))))))
    svg = atlasPage.drawQuiver(data)
    assert svg.startswith('<svg') and svg.count('<circle') == 5 and svg.count('<path') >= 4


def test_the_page_carries_every_section(atlas5):
    page = atlasPage.render('Shape atlas n = 5', atlasPage.sectionsFrom(atlas5['tables'], 2, 3))
    assert page.startswith('<!doctype html>') and '<title>Shape atlas n = 5</title>' in page
    assert page.count('<svg') >= 3 and 'prefers-color-scheme: dark' in page


def test_the_command_line_reads_a_ledger(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    args = argparse.Namespace(length = 5, depth = 2, sample = 0, seed = 0)
    jobs.runTask(batch.TASKS['atlas'], args, out = io.StringIO())
    assert atlas.main(['5', '--depth', '2', '--validate', '--page', 'page.html', '--top', '3']) == 0
    printed = capsys.readouterr().out
    assert 'hubs' in printed and 'candidate merges: 0' in printed and 'validation' in printed
    assert (tmp_path / 'page.html').exists()
    assert (tmp_path / 'logs' / 'atlas-n5-d2-s0-r0.nodes.parquet').exists()
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `$WSL .venv/bin/python -m pytest tests/test_shape_atlas.py -q`
Expected: FAIL at collection with `ModuleNotFoundError: No module named 'atlas'`.

- [ ] **Step 3: Write `quivermutation/atlasPage.py`**

```python
"""The atlas's top shapes as a page to look at, not keys to read.

`classpage` draws lines and quipus in the browser from their names; the shapes
here have no names, so they are drawn in Python from the quivers themselves --
a spring layout with a fixed seed, parallel arrows bowed apart, and each
relation written out underneath.  The page is one static file with everything
inline and follows the light and dark schemes of whatever opens it.
"""

import html
import math

import networkx as nx
import polars as pl

from . import shapeAtlas
from . import shapeKeys


_STYLE = """
:root { --bg: #fbfaf7; --fg: #1d1d1b; --muted: #6b6a66; --line: #3b3a36;
        --node: #ffffff; --accent: #b4532a; --card: #f2f0ea; }
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) { --bg: #161615; --fg: #ecebe6; --muted: #9c9a93;
        --line: #cfcdc5; --node: #262624; --accent: #e08a5f; --card: #1f1f1d; }
}
:root[data-theme="dark"] { --bg: #161615; --fg: #ecebe6; --muted: #9c9a93;
        --line: #cfcdc5; --node: #262624; --accent: #e08a5f; --card: #1f1f1d; }
body { background: var(--bg); color: var(--fg); margin: 0 auto; max-width: 1100px;
       padding: 24px 16px; font: 15px/1.5 system-ui, sans-serif; }
h1 { font-size: 1.5rem; } h2 { font-size: 1.15rem; margin-top: 2rem; }
.grid { display: grid; gap: 16px; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); }
.card { background: var(--card); border-radius: 8px; padding: 12px; overflow-wrap: anywhere; }
.card h3 { font-size: .95rem; margin: 0 0 6px; }
.card p { color: var(--muted); font-size: .8rem; margin: 2px 0; }
svg { width: 100%; height: auto; }
"""


def drawQuiver(data, size = 240):
    """One quiver as an inline SVG: vertices labelled, arrows with heads."""
    algebra = shapeKeys.deserialise(data)
    quiver = algebra.quiver
    graph = nx.Graph(quiver.to_undirected())
    if graph.number_of_nodes() > 1:
        raw = nx.spring_layout(graph, seed = 0)
    else:
        raw = {vertex: (0.0, 0.0) for vertex in graph.nodes}
    xs = [p[0] for p in raw.values()] or [0.0]
    ys = [p[1] for p in raw.values()] or [0.0]
    margin = 22
    span = max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
    scale = (size - 2 * margin) / span
    at = {v: (margin + (x - min(xs)) * scale, margin + (y - min(ys)) * scale)
          for v, (x, y) in raw.items()}
    parts = ['<svg viewBox="0 0 {0} {0}" xmlns="http://www.w3.org/2000/svg" role="img">'.format(size),
             '<defs><marker id="head" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
             'markerHeight="6" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="var(--line)"/>'
             '</marker></defs>']
    bundles = {}
    for tail, head, key in sorted(quiver.edges(keys = True)):
        bundles.setdefault((tail, head), []).append(key)
    radius = 9
    for (tail, head), keys in bundles.items():
        (x1, y1), (x2, y2) = at[tail], at[head]
        length = math.hypot(x2 - x1, y2 - y1) or 1.0
        ux, uy = (x2 - x1) / length, (y2 - y1) / length
        sx, sy = x1 + ux * radius, y1 + uy * radius
        ex, ey = x2 - ux * (radius + 2), y2 - uy * (radius + 2)
        for position, _key in enumerate(keys):
            bend = (position - (len(keys) - 1) / 2) * 18
            cx, cy = (sx + ex) / 2 - uy * bend, (sy + ey) / 2 + ux * bend
            parts.append('<path d="M{0:.1f},{1:.1f} Q{2:.1f},{3:.1f} {4:.1f},{5:.1f}" '
                         'fill="none" stroke="var(--line)" stroke-width="1.4" '
                         'marker-end="url(#head)"/>'.format(sx, sy, cx, cy, ex, ey))
    for vertex, (x, y) in at.items():
        parts.append('<circle cx="{0:.1f}" cy="{1:.1f}" r="{2}" fill="var(--node)" '
                     'stroke="var(--accent)" stroke-width="1.4"/>'.format(x, y, radius))
        parts.append('<text x="{0:.1f}" y="{1:.1f}" font-size="9" text-anchor="middle" '
                     'fill="var(--fg)">{2}</text>'.format(x, y + 3, html.escape(str(vertex))))
    parts.append('</svg>')
    return ''.join(parts)


def render(title, sections):
    """A whole page: a heading per section, a card per shape."""
    body = ['<h1>{0}</h1>'.format(html.escape(title))]
    for heading, entries in sections:
        body.append('<h2>{0}</h2><div class="grid">'.format(html.escape(heading)))
        for entry in entries:
            body.append('<div class="card"><h3>{0}</h3>{1}{2}</div>'.format(
                html.escape(entry['heading']), drawQuiver(entry['quiver']),
                ''.join('<p>{0}</p>'.format(html.escape(line)) for line in entry['lines'])))
        body.append('</div>')
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1">'
            '<title>{0}</title><style>{1}</style></head><body>{2}</body></html>'
            .format(html.escape(title), _STYLE, ''.join(body)))


def sectionsFrom(tables, level, top):
    """Hubs, bridges and leftover hubs at one level, as page sections."""
    import json
    key = 'key{0}'.format(level)
    measures = shapeAtlas.shapeMeasures(tables, level).filter(~pl.col('isLine'))
    nodes = tables['nodes']

    def entries(frame):
        found = []
        for row in frame.head(top).iter_rows(named = True):
            quiver = json.loads(nodes.filter(pl.col(key) == row[key])['quiver'][0])
            found.append({
                'heading': '{0} classes, {1} starts'.format(row['classes'], row['starts']),
                'quiver': quiver,
                'lines': [shapeKeys.describe(shapeKeys.deserialise(quiver)),
                          'return rate {0:.2f}, first at depth {1}, leftover share {2:.2f}'.format(
                              row['returnRate'], row['firstDepth'], row['leftoverShare'])],
            })
        return found

    return [
        ('Hubs', entries(measures)),
        ('Bridges', entries(measures.filter(pl.col('bridge')))),
        ('Hubs among leftovers', entries(measures.filter(pl.col('leftoverShare') > 0)
                                         .sort('leftoverShare', 'starts', descending = True))),
    ]
```

- [ ] **Step 4: Write `atlas.py`**

```python
#!/usr/bin/env python
"""Read a shape-atlas ledger: which quiver shapes the walks pass through.

    python batch.py atlas 9 --depth 4 --jobs 7          walk, and write the ledger
    python atlas.py 9 --depth 4                         read it
    python atlas.py 9 --depth 4 --validate --page logs/atlas-n9.html

Resolves every quiver in the ledger into exact keys at four levels up to
relabelling (`quivermutation/shapeKeys.py`), writes the tables as parquet beside
the ledger, and prints the hubs, the bridges, the hubs among the leftovers, the
commonest cycles through a line, and every candidate merge -- a shape two
classes share -- after replaying it.  `--validate` runs H-022's four checks.
Spec: docs/superpowers/specs/2026-09-24-shape-atlas-design.md.
"""

import argparse
import sys

from quivermutation import atlasPage
from quivermutation import jobs
from quivermutation import shapeAtlas


def main(argv = None):
    parser = argparse.ArgumentParser(description = __doc__,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("length", type = int)
    parser.add_argument("--depth", type = int, default = 4)
    parser.add_argument("--sample", type = int, default = 0)
    parser.add_argument("--seed", type = int, default = 0)
    parser.add_argument("--level", type = int, default = 2, choices = (0, 1, 2, 3),
                        help = "the level hubs and cycles are counted at (default 2)")
    parser.add_argument("--top", type = int, default = 30)
    parser.add_argument("--validate", action = "store_true",
                        help = "run H-022's four checks")
    parser.add_argument("--page", default = None, help = "write the drawings here")
    args = parser.parse_args(argv)

    path = shapeAtlas.ledgerPath(args.length, args.depth, args.sample, args.seed)
    records = jobs.Ledger(path).records()
    if not records:
        print("no ledger at {0}; run `python batch.py atlas {1} --depth {2}` first".format(
            path, args.length, args.depth))
        return 1
    tables = shapeAtlas.resolve(records, args.length)
    shapeAtlas.writeTables(tables, path[:-len('.jsonl')])
    shapeAtlas.report(tables, args.level, args.top, sys.stdout)

    candidates = shapeAtlas.candidateMerges(tables)
    print("\ncandidate merges: {0}".format(len(candidates)))
    for candidate in candidates:
        verdict = shapeAtlas.replay(args.length, candidate)
        print("  {0} {1} ~ {2}  via {3} / {4}: {5}".format(
            'MERGE ' if verdict['ok'] else 'FAILED', candidate['first']['cls'],
            candidate['second']['cls'], candidate['first']['path'],
            candidate['second']['path'], verdict['reason']))

    if args.validate:
        print("\nvalidation (H-022)")
        for name, outcome in shapeAtlas.validate(tables, args.length, records).items():
            print("  {0}: {1}".format(name, outcome))
    if args.page:
        title = "Shape atlas n = {0}, depth {1}".format(args.length, args.depth)
        with open(args.page, "w", encoding = "utf-8") as handle:
            handle.write(atlasPage.render(title, atlasPage.sectionsFrom(
                tables, args.level, args.top)))
        print("\nwrote {0}".format(args.page))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 5: Register the script**

In `pyproject.toml`, change the `py-modules` line to:

```toml
py-modules = ["atlas", "batch", "classify", "classes", "discover", "families", "merges", "overlaps", "overnight", "probe"]
```

In `batch.py`, append to the `ELSEWHERE` list (after the `overlaps` entry):

```python
    ("atlas", "python atlas.py 9 --depth 4 --validate --page logs/atlas-n9.html",
     "read a `batch.py atlas` ledger: hubs, bridges, cycles, candidate merges"),
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `$WSL .venv/bin/python -m pytest tests/test_shape_atlas.py tests/test_shape_keys.py -q`
Expected: all pass.

- [ ] **Step 7: Run the whole fast suite**

Run: `$WSL .venv/bin/python -m pytest -q -m "not slow"`
Expected: all pass (nothing outside the new files changed behaviour; `batch.py`'s only change is the new task and one `ELSEWHERE` line).

- [ ] **Step 8: Commit**

```bash
git add quivermutation/atlasPage.py atlas.py pyproject.toml batch.py tests/test_shape_atlas.py
git commit -m "Read an atlas ledger from the command line, and draw its top shapes

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Write H-022, the glossary terms and the README section — before any run

**Files:**
- Modify: `research/HYPOTHESES.md` (new entry at the top, directly after the header's `---`)
- Modify: `GLOSSARY.md` (a new `## Shapes` section before the `## Equivalences that save work` section; add it to the Contents list)
- Modify: `README.md` (a subsection at the end of "Other families that could carry the classes the theorem misses")

- [ ] **Step 1: Add H-022 to `research/HYPOTHESES.md`**

Insert directly above `## H-021 — ...`:

```markdown
## H-022 — The walks pass through a small set of shapes, and the ones two classes share are merges
*2026-09-24, written before the first census* · **OPEN**

Three non-line families have mattered so far -- quipus with relations (F-034),
the squares with a side of two (F-027), one parallel pair (H-016) -- and each
was found by reading one walk by hand. The **shape atlas** (`batch.py atlas`,
`atlas.py`; spec `docs/superpowers/specs/2026-09-24-shape-atlas-design.md`)
records every quiver the depth-4 walks reach out of every LNA of `n = 8` and
`9`, and out of the `n = 10` leftovers with 20 of each quipu class, keyed up to
relabelling at four levels (L0 graph, L1 quiver, L2 relation skeleton, L3
algebra).

**The hypothesis.** A small number of shapes carries most of the walks between
lines -- hubs -- and a shape reached from two classes is a merge the label-exact
meeting of `search.meetingPoints` cannot see, because it compares quivers with
their labels.

**What must come out first, or the instrument is wrong** (`atlas.py --validate`):

1. Among the non-line L2 shapes that lead back to a line, the commonest with a
   square has a short side of **two**, and no square has a short side of three
   (F-027).
2. At `n = 9`, at least seven of the nine leftover LNAs reach `P^(6)_(1,1)`
   with relations (H-014).
3. On every quipu with monomial relations, equal L3 keys and equal
   `quipuRelations.certificate` coincide exactly.
4. Everything the label-exact search reaches at depth 2 is in the census.

A failure of 1 or 2 means nothing else the atlas says is read. A failure of 3
or 4 is a bug in the keys.

**What would settle it.** Yes: a replayed candidate merge between two orbits
that no move and no label-exact search has joined -- at `n = 10`, between two
of H-013's leftover orbits. No: every shape shared across classes is shared
only by classes already known to be one, at depth 4 and at depth 5, which says
the non-line shapes are a detour at these lengths and not a shortcut. Either is
worth having.

---
```

- [ ] **Step 2: Add the terms to `GLOSSARY.md`**

Add `- [Shapes](#shapes)` to the Contents list in the position matching where the section goes, and insert this section directly above `## Equivalences that save work`:

```markdown
## Shapes

**Shape.** A quiver with relations up to relabelling its vertices, at one of
four levels of detail (`quivermutation/shapeKeys.py`):

| level | keeps |
|---|---|
| **L0** | the underlying undirected multigraph |
| **L1** | the quiver, parallel arrows included |
| **L2** | the quiver and its **relation skeleton**: each relation's kind (zero, commutativity, other) and path lengths, on the vertices it runs through |
| **L3** | the algebra as presented, up to relabelling, parallel-arrow naming and the sign gauge |

**Bucket.** The Weisfeiler–Lehman hash of a level's graph. Isomorphic shapes
share a bucket; so, rarely, do non-isomorphic ones, which is why a bucket is
resolved into **keys** by an exact isomorphism test (`ShapeIndex`).

**Label-exact.** Equal with the labels as they stand: `search.quiverKey`,
`fingerprint.canonicalKey`. Right inside one walk, where labels do not move
(F-049); blind to two walks reaching one quiver under different labels.

**Shape atlas.** The census of every quiver the walks out of a set of LNAs
reach (`batch.py atlas`), read by `atlas.py`. H-022.

**Hub.** A shape reached from many classes. **Bridge.** A shape reached from
two or more orbits the classification already puts in one class.
**Return rate.** Of the starts reaching a shape, the share whose walk goes on
from it to a line other than the start.

**Candidate merge.** An L3 shape reached from two classes. It is a merge only
once `shapeAtlas.replay` has re-run both paths under the Coxeter guard and
found the ends isomorphic.
```

- [ ] **Step 3: Add the README section**

In `README.md`, directly above the `### Reorienting a tree is free` heading, insert:

```markdown
### Which shapes the walks pass through

```bash
python batch.py atlas 9 --depth 4 --jobs 7            # walk every LNA, record every quiver
python atlas.py 9 --depth 4 --validate --page logs/atlas-n9.html
python batch.py atlas 10 --depth 4 --sample 20 --jobs 7   # the leftovers and 20 per quipu class
```

Every family above was found by reading one walk by hand. The atlas counts
instead: it records every quiver the walks reach, keyed up to relabelling at
four levels -- graph, quiver, relation skeleton, algebra -- and reports the
**hubs** (shapes many classes pass through), the **bridges** (shapes where the
classification already joined two orbits), the commonest cycles line → shape →
line, and every shape two classes share, replayed before it is called a merge.
Label-exact meeting cannot see those: two walks out of different LNAs reach the
same quiver under different labels. See research H-022 and `GLOSSARY.md`,
"Shapes".
```

- [ ] **Step 4: Commit**

```bash
git add research/HYPOTHESES.md GLOSSARY.md README.md
git commit -m "Write down what the shape atlas must find before it is believed (H-022)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Run the census at n = 8, 9 and 10, and record it

**Files:**
- Modify: `research/EXPERIMENTS.md` (E-053 at the top), `research/HYPOTHESES.md` (H-022 status line)
- Possibly create: `research/FINDINGS.md` F-054, if the runs establish something

This task is a run, not code. Do each length in order and **stop at the first failed validation check 1 or 2** — the instrument is then wrong, and that goes back to Task 4 or 5 as a bug, not into the research record as a finding.

- [ ] **Step 1: n = 8, every LNA**

```bash
MSYS_NO_PATHCONV=1 wsl.exe -d Ubuntu --cd /mnt/c/Users/didri/kode/quiverMutation -- .venv/bin/python batch.py atlas 8 --depth 4 --jobs 7
MSYS_NO_PATHCONV=1 wsl.exe -d Ubuntu --cd /mnt/c/Users/didri/kode/quiverMutation -- .venv/bin/python atlas.py 8 --depth 4 --validate --page logs/atlas-n8.html > logs/atlas-n8.txt
```

Expected: 429 starts, `candidate merges: 0` (n = 8 is fully classified, so any candidate is a key bug or a search bug — investigate before going on), `certificates` and `coverage` ok, `quipuHub: skipped`. Record the wall clock of both commands.

- [ ] **Step 2: n = 9, every LNA**

Same two commands with `9` and `logs/atlas-n9.*`. Expected: 1430 starts; the validation `squares` and `quipuHub` checks ok. If `quipuHub.leftoverStarts` is 0, first check the name: print `features['quipu']` for the tree nodes the leftover `3033030` reaches and compare with `N9_HUB` — H-014's name came from `families.py members`, and a naming mismatch is not a refutation.

- [ ] **Step 3: n = 10, leftovers and 20 per quipu class**

```bash
MSYS_NO_PATHCONV=1 wsl.exe -d Ubuntu --cd /mnt/c/Users/didri/kode/quiverMutation -- .venv/bin/python batch.py atlas 10 --depth 4 --sample 20 --jobs 7
MSYS_NO_PATHCONV=1 wsl.exe -d Ubuntu --cd /mnt/c/Users/didri/kode/quiverMutation -- .venv/bin/python atlas.py 10 --depth 4 --sample 20 --validate --page logs/atlas-n10.html > logs/atlas-n10.txt
```

Look in particular at the candidate merges between leftover orbits (`orbit:` classes): H-013's groups at n = 10 are where a merge would be new. Every `MERGE` line is a claim; every `FAILED` line is a bug to explain.

- [ ] **Step 4: Record E-053**

Insert at the top of `research/EXPERIMENTS.md`, above `## E-052`, filling each bracket from the three runs' output files (`logs/atlas-n*.txt`) — every number comes from a run, none is estimated:

```markdown
## E-053 — The first shape atlas: n = 8 and 9 whole, n = 10 leftovers, depth 4
*2026-09-24* · tests H-022

**Commands.** `python batch.py atlas {8,9} --depth 4 --jobs 7`,
`python batch.py atlas 10 --depth 4 --sample 20 --jobs 7`, each read with
`python atlas.py <n> --depth 4 [--sample 20] --validate --page logs/atlas-n<n>.html`.

| n | starts | quivers recorded | distinct (labels) | L0 / L1 / L2 / L3 shapes | census wall clock | analysis wall clock |
|---|---|---|---|---|---|---|
| 8 | [..] | [..] | [..] | [..] | [..] | [..] |
| 9 | [..] | [..] | [..] | [..] | [..] | [..] |
| 10 | [..] | [..] | [..] | [..] | [..] | [..] |

**Validation (H-022).** squares: [short sides counted, top square]; quipu hub:
[leftover starts reaching P^(6)_(1,1)]; certificates: [checked, violations];
coverage: [starts checked, missing].

**Top hubs at L2**, per length: the first five from each report, with their
`describe` line and classes / starts / return rate.

**Candidate merges.** [count per length; for each MERGE the two classes and the
two paths; for each FAILED the reason].

**What it cost and what to run next.** [one paragraph].
```

- [ ] **Step 5: Update H-022's status line**

Change H-022's status from `**OPEN**` to what the runs support. That is `**SUPPORTED**` if the four checks pass and hubs are concentrated (say what "concentrated" was measured as). If a replayed merge joined two orbits nothing had joined before, add a finding `F-054` in `research/FINDINGS.md` with the pair, both paths and the replay, and point H-022 at it (`**CONFIRMED → F-054**` only if that is the whole claim). Put the evidence in brackets after the status, as the other entries do, citing E-053.

- [ ] **Step 6: Commit**

```bash
git add research/EXPERIMENTS.md research/HYPOTHESES.md research/FINDINGS.md
git commit -m "Run the first shape atlas at n = 8 to 10, and record what it found (E-053)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
