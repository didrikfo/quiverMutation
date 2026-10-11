# Shape atlas — design

*2026-09-24*

## Why

The classification has been driven by walking mutations out of LNAs and merging
what meets. Three families of non-line quivers have turned out to matter along the
way — quipus with relations (F-034), the commutative squares with a side of two
(F-027), and quivers with one parallel pair (H-016) — and each was found by
looking at one walk by hand. Nothing yet asks the question the other way round:
**over every walk, which quiver shapes do the walks pass through, and which of
those tie classes together?**

There is also a concrete gap that makes the question impossible to ask today.
`search.quiverKey` is exact *with labels*, and returns `None` for any quiver with
parallel arrows. Inside one walk that is correct — labels do not move (F-049) —
but two walks out of different LNAs that reach **isomorphic, differently
labelled** quivers are never seen to meet. The atlas needs a relabelling-
invariant key, and once it has one, every shape shared between two classes is
also a candidate merge.

This is a survey. Its output is data and a first reading of it, from which the
next line of work (hub landmarks, a shape-level move grammar, theory-supplied
families such as canonical algebras) is chosen.

## What it is

A census task and an analysis over its ledger.

1. **Census** — walk to a fixed depth out of a chosen set of starting LNAs,
   record every quiver reached under four keys, and every mutation step between
   them.
2. **Analysis** — per shape, how widely it is reached, how much it mixes
   classes, how early it appears, and whether walks through it come back to a
   line; the transition structure between shapes; candidate merges at the finest
   level, verified before they are reported.
3. **Page** — the top hubs and bridges drawn, so shapes can be looked at rather
   than read as keys.

## The four keys

Every quiver reached gets four keys, coarse to fine. All are invariant under
relabelling the vertices.

| level | what it forgets | what it keeps |
|---|---|---|
| **L0** | orientation, relations | underlying undirected multigraph |
| **L1** | relations | the quiver, parallel arrows included |
| **L2** | coefficients, which terms | the quiver plus a *relation skeleton*: each relation as its kind (`zero` for a monomial, `comm` for a two-term relation, `other` for anything longer) and the multiset of its path lengths, attached as a hyperedge on the vertices its paths run through, with its start and end marked |
| **L3** | nothing but labels | the algebra: quiver plus the relations as presented, with parallel-arrow naming and the sign gauge quotiented as `fingerprint.canonicalKey` does |

**Computing them.** Each level builds a labelled graph — vertices, one node per
arrow, and for L2/L3 one node per relation and per path term, with edge labels
carrying order along a path — and takes
`networkx.weisfeiler_lehman_graph_hash` of it as a **bucket**. A bucket is not a
key: WL can put non-isomorphic graphs in one bucket. Within a bucket, quivers are
compared exactly:

- L0–L2: an exact isomorphism test (`networkx` VF2 matchers with node and edge
  labels) between the new quiver and each representative already in the bucket.
  The key is `bucket:index-of-representative`.
- L3: for each isomorphism the L2 matcher finds, relabel the algebra through it
  and compare `fingerprint.canonicalKey`. Equal on some isomorphism means the
  same algebra up to relabelling. This is **sound but conservative** in the way
  F-049 already records: two presentations of one ideal by different generators
  are two keys. A false negative loses a meeting; there are no false positives.

Representatives are held per bucket in the analysis, not the census: the census
writes the WL bucket plus a serialised quiver for each distinct node, and the
analysis resolves buckets into exact keys once, over the whole ledger. That keeps
`run` a pure function of its start (a `jobs.Task` requirement) and keeps every
exact comparison in one place.

The existing `quipuForms.canonicalTreeForm` and the quipu certificate of
`quipuRelations` are exact on trees and are reused there as a cross-check: on
every tree node, equal L3 keys must mean equal certificates and vice versa.

## Census

**Task.** `batch.py atlas <n>`, a `jobs.Task` like `sample` and `cores`:
append-only ledger under `logs/`, one line per finished start, resumable,
`--jobs`, `--budget-hours`, exit 2 on budget. The ledger name carries `n`, the
depth and the start set, since those change the work; `--jobs` and filters do
not.

**Starts.**

| n | starts | why |
|---|---|---|
| 8 | every LNA (429) | complete, cheap, and fully classified (F-011), so every class label is certain |
| 9 | every LNA (1430) | complete, and holds the first two classes outside the theorem |
| 10 | every LNA outside a quipu class (262), plus up to 20 per quipu class drawn uniformly with a fixed seed | the leftovers are the point; the sample gives the contrast |

Each start is walked **and** its relation dual, as `quiversReachedFrom` does,
with what the dual reaches carried back through the opposite so every key is a
quiver the start itself reaches.

**Tag per start.** Its orbit under `freeMoves.derivedOrbits(n, free = True,
edges = True, doubles = True)` closed under the dual, and its class where the
classification names one (quipu name, canonical type, or the orbit key for a
leftover). Class is what hub scores count over; orbit is kept so a shape that
joins two orbits of one class can be told apart from one that joins two classes.

**Depth.** 4 for the first run. Measured on 2026-09-24, one direction, with
`fingerprint.Visited`: `3345000` 126 distinct nodes in 1.6 s, `3033030` 110 in
1.8 s, `34504030` (n = 10) 138 in 3.0 s, `0000000` 558 in 2.7 s. Both directions
at ~5 s a start gives about 40 min for n = 8 and 9 together and 30 min for n = 10
on one core, so the whole census is well under an hour on seven. Depth 5 is the
second run, on whatever depth 4 shows to be worth it.

**What a ledger line holds.**

```
start       relation lengths, e.g. "3345000"
orbit, cls  the tags above
nodes       [ {bucket0, bucket1, bucket2, bucket3, depth, features, quiver} ... ]
edges       [ [parentIndex, vertex, childIndex] ... ]    (vertex signed: + right, - left)
```

`quiver` is the node serialised (arrows with keys, relations as `arrowRels`) so
the analysis can rebuild it. `depth` is the shortest depth the node was reached
at from this start. Edges come from the visitor's path: the parent of a node
reached by `path` is the node recorded at `path[:-1]`, which the depth-first walk
visits first. `features` are cheap and computed in the worker:

- vertex, arrow and relation counts; number of parallel bundles; undirected
  cycle rank;
- `is_line`, `is_tree`, `is_quipu` (by `quipuForms.isQuipuByDegrees` on trees);
- `defect` = relations − cords (H-017), on quipus only;
- Coxeter polynomial, as a string.

## Analysis

`atlas.py <n>` reads a ledger, resolves exact keys, writes three parquet tables
beside it (`nodes`, `visits`, `edges`, one row per distinct key / per
start-and-key / per distinct key-to-key step at each level) and prints a report.
`--page FILE` writes the drawings.

**Per shape, at each of L0–L3:**

| measure | meaning |
|---|---|
| reach | number of starts, orbits and classes whose walk passes through it |
| mixing | number of classes reached through it, and the entropy of the start distribution over classes |
| first depth | median and minimum depth it is first reached at |
| return rate | of the starts reaching it, the fraction whose walk reaches a line *below* it within the remaining depth |
| leftover share | fraction of its starts that are outside a quipu class |

**Hubs** are the shapes ranked by class reach. **Bridges** are shapes reached
from two orbits that the classification says are one class — these are the
shapes that tied the classes together. At L2, the report lists the top 30 of
each and, separately, the top shapes among leftovers only.

**Shape transitions.** Collapse `edges` to L2 keys and count. Report the most
common line → S → … → line cycles of length up to 4 in the collapsed graph:
these are rule *templates* in the sense of F-027, and the squares should be most
of the length-2 ones.

**Candidate merges.** An L3 key reached from starts in two different orbits —
orbits that are *not* already one class — is a candidate merge, carried with a
path from each side. Before it is printed it is replayed: both paths re-run step
by step under the Coxeter guard, the endpoint rebuilt, the isomorphism
re-derived, and the two starts' Coxeter polynomials compared. Only replayed
candidates are reported as merges; the rest are listed as failed replays, which
would be a bug in the key and must be looked at, not dropped.

## Page

A static page like `classes.py --page`: for each of the top hubs and bridges at
L2, one representative drawn as a quiver with its relations, its measures, and
the classes it touches. Reuses the drawing code in `classpage`, extended to
quivers that are not lines (layout by `networkx` spring or planar layout; the
existing code draws only lines and quipus).

## Validation — before trusting any of it

Written into the research harness as **H-022** before the run, with these as its
predictions:

1. **The squares.** Among L2 shapes on walks that return to a line within the
   depth, the commutative squares with a short side of two are the most common
   non-line shape, and no square with a short side of three appears (F-027).
2. **The n = 9 quipu.** `P^(6)_(1,1)` with relations is a top hub at L1 for the
   nine n = 9 leftovers, reached by at least the seven H-014 names.
3. **Keys agree with what is known.** On every tree node, L3 keys and quipu
   certificates agree exactly. On every start, the lines its walk reaches under
   the old label-exact key are a subset of those reached under L3.
4. **Known meetings come back.** For every pair of n = 8 starts that
   `search.meetingPoints(first, second, 2)` joins under the label-exact key, the
   census at depth 4 has a shared L3 key between the two, since the label-exact
   meeting is one.

A failure of 1 or 2 means the instrument is wrong and nothing else it says is
read. A failure of 3 or 4 is a bug in the keys.

## Tests

- **Keys** (`tests/test_shape_keys.py`): relabelling a quiver with relations by
  a random permutation leaves all four keys unchanged; two non-isomorphic
  quivers that WL hashes together (a small constructed pair) get different
  exact keys; a commutativity relation and a zero relation on the same square
  differ at L2; the two namings of a parallel bundle and a sign flip of an arrow
  agree at L3.
- **Census** (`tests/test_atlas_task.py`): the task at n = 5, depth 2, on three
  starts — ledger lines round-trip, resume skips done starts, a depth-2 walk's
  node set equals the one `quiversReachedFrom` gives under the old key wherever
  that key is defined.
- **Analysis**: at n = 5, every LNA is one of a handful of known classes; the
  report's class reach for the lines themselves matches the classification.
  Candidate-merge replay is exercised on one planted pair.

All fast. The n = 8 to 10 census is a run, not a test.

## Research harness

- **H-022** — the atlas hypothesis and the four predictions above, written
  before the run.
- **E-0nn** — each census and analysis run, with the command and what it cost.
- **F-0nn** — whatever it finds, including a clean negative: "no non-line shape
  is shared across classes at depth 4" is worth having.
- **GLOSSARY.md** — *shape*, *L0–L3*, *hub*, *bridge*, *return rate*.

## Out of scope

Chosen later, from what the atlas shows:

- a landmark index — deep walks out of the top hubs, joined once and reused by
  every LNA;
- a shape-level move grammar and best-first search or pruning built from it;
- theory-supplied families (canonical and extended canonical algebras, the
  vect-X(2,a,b) route of F-046) as seeds for the census;
- n = 11 and longer, and depth 5 or more.

The label-exact `quiverKey` stays as it is. Replacing it inside `merges.py` is
a follow-up once the L3 key has been validated here.
