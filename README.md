# quiverMutation

Tools for performing mutations of quivers with relations, corresponding to tilting
mutations of their path algebras.

The mutation procedure implemented here is the combinatorial rule for tilting
mutation from:

* D. Fosse, *A combinatorial procedure for tilting mutation*,
  [arXiv:2112.08129](https://arxiv.org/abs/2112.08129)

The main application is the classification of linearly oriented Nakayama algebras
(LNAs) up to derived equivalence, which produced the results in:

* D. Fosse, *Quipu quivers and Nakayama algebras with almost separate relations*,
  [arXiv:2305.06642](https://arxiv.org/abs/2305.06642)

## Getting started

With [uv](https://docs.astral.sh/uv/):

```bash
uv venv --python 3.11
uv pip install -e '.[test]'
```

Or with plain venv/pip:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[test]'
```

Run the test suite:

```bash
.venv/bin/python -m pytest -q
```

## Where things are written down

* [`NOTES.md`](NOTES.md) — the code: what the model expresses, known gaps, the backlog.
  The package layout is listed in [`quivermutation/__init__.py`](quivermutation/__init__.py).
* [`research/`](research/) — the mathematics: findings, hypotheses, retractions,
  the log of runs made, and summaries of the literature. All dated, nothing
  deleted. Read [`research/README.md`](research/README.md) before adding to it.

## Classifying a length

```bash
python classify.py 8
```

This classifies all 429 LNAs on 8 vertices up to derived equivalence and writes
`A_8_mutation_classes.csv` (and a parquet alongside it) into the current
directory, with one row per LNA giving

| column | meaning |
|---|---|
| `Relations` | the LNA, as `1;2;3\|3;4;5;6` -- one `;`-joined path per relation |
| `Mutation class` | the class, named by its quipu |
| `Mutation path from class representative` | the mutations that get there |
| `Coxeter polynomial` | a derived invariant, though not a complete one |
| `Numbering` | the vertex numbering the mutation path produces |
| `Hereditary form` | the quipu the class corresponds to |

It prints the classes and their sizes, and exits non-zero if any class was left
unsettled.

A long run does not have to finish in one sitting. Every step writes the table
after every class and records what it finished in a JSON file beside the CSV, so
`--resume` picks up where it stopped rather than redoing the naming and the
resolving:

```bash
python classify.py 10 --resume
```

`--budget-hours` stops a run cleanly once the budget is spent, between classes,
with everything done so far on disk -- which is how to fit a classification into
a fixed window such as a night. It exits 2 when it stops that way, so a wrapper
can tell "out of time, resume me" from "finished, with something unsettled":

```bash
python classify.py 10 --budget-hours 9 --resume
```

[`overnight.sh`](overnight.sh) is that wrapper: it keeps the machine awake,
restarts a job that dies, and runs the deep interior probe of research H-010
alongside the classification.

```bash
./overnight.sh 9
```

The classification runs in four steps, described in `NOTES.md`: seed every LNA
the quipu theorem of arXiv:2305.06642 covers, search by mutation for the rest,
name any class the theorem missed by the hereditary algebra its search reaches,
and settle whatever is left by a deeper search. For n <= 8 this reproduces the
published classification with nothing left over, replacing what used to be a
hand-merge over the CSV.

## Reading a classification back

A table with one row per LNA is the wrong shape for looking at the answer, and
there are Catalan(n-1) of them -- 1430 at n = 9, 58786 at n = 12 -- while the
number of classes stays small. `classes.py` reads the table by class instead:

```bash
python classes.py 9                       # the classes, largest first
python classes.py 9 --collisions          # where the Coxeter polynomial stops separating
python classes.py 9 --kind "not piecewise hereditary"
python classes.py 9 --members "P^(1,4)_(1,0,1)"
python classes.py 9 --page A_9.html       # the same thing as a page to browse
```

Everything but `--members` reads the columns it needs and groups; nothing loads
the per-LNA rows for the whole table. The page carries the classification inline
-- no server, nothing to fetch -- and draws each class' quipu and each LNA as its
quiver with an arc over the span of every relation.

The same thing from Python:

```python
from quivermutation import classview

nine = classview.Classification.forLength(9)
nine.classes(kind = classview.QUIPU, minSize = 100)   # the big quipu classes
nine.coxeterCollisions()                              # what the polynomial cannot separate
nine.members("P^(1,4)_(1,0,1)")                       # one class, with the path to each member
```

## Where a classification still has to search

```bash
python overlaps.py 6 7 8 9              # coverage by relation overlap
python overlaps.py 9 --cores            # what is left, by overlapping run
python overlaps.py 9 --free             # with relations of two arrows free
python probe.py 1:3,2:3 --steps 4       # what one configuration can become
```

The quipu theorem names the class of an LNA whose consecutive relations share at
most one arrow, and the move rules of `lnaMoves` carry the rest into its reach --
100% of them at n = 6 and n = 7, 98% at n = 8, 84% at n = 9, 63% at n = 10 and
47% at n = 11, with no search run at all. What is left over is exactly the LNAs
with two relations sharing two or more arrows, and `overlaps.py` prints where
the boundary sits at each length, which configurations are stuck, and how much
each half of the rule table is worth. See `research/` F-021 to F-025.

`probe.py` is the other half of the same question: instead of "what rules are
there", it asks what can happen to one named configuration, and reports what it
reaches grouped by overlap. A run that reaches nothing with a smaller overlap is
a negative result worth having -- that is how the obstruction above was found.

The two halves are different in kind. A *floating* rule holds at every position
of the quiver; an *anchored* one holds only against the source or the sink, where
a mutation does something it cannot do in the interior. The anchored half does
most of the work above the almost separate line. Both are closed under the
relation dual -- reverse every arrow and exchange right mutation for left, which
takes a rule to a rule (`lnaMoves.dualRule`).

Two things reach further than any rule does, and neither is a table row.
`freeMoves` deletes a relation of **two arrows**, which arXiv:2310.08346 says
leaves the derived equivalence class alone -- no mutation, no sequence, and it
merges more at n = 12 than the whole rule table does. `edgeMoves` holds a family the rule
encoding cannot state at all: a relation at an end of the quiver doubles, and its
window is allowed to be crossed by a relation it never touches. With both,
**n = 8 needs no search at all** -- 21 orbits, nothing left over -- and n = 9
falls from 222 rows to 37. See `research/` F-028 to F-030.

Past `n = 12` the orbits cannot be partitioned at all -- `n = 13` has 208012 LNAs
and the interesting rows are a handful -- so `freeMoves.orbitOf` walks the moves
out of one row, and `freeMoves.movesJoin` walks out of **two rows at once** and
stops where they meet. Use the second for a membership question: a one-way walk
that stops at its row cap has measured the budget and not the moves, which is how
49 barricades at `n = 15` and `16` were first recorded as failures and then joined
in 45 seconds. See `research/` F-042 and E-037.

## Lengths too long to enumerate

Past about `n = 13` there is no complete pass to be had -- 208012 LNAs at
`n = 13`, 1767263190 at `n = 20` -- and the lengths that *can* be done completely
are unrepresentative, because in a quiver of length 8 every vertex is within
three arrows of an end and the anchored rules are exactly the ones an end is in
reach of. So the question a long length can be asked is a statistical one:

```bash
python batch.py sample 16 --count 4000 --jobs 7 --budget-hours 9
python batch.py sample 16 --summary
```

This draws LNAs **uniformly** -- exactly so, by a dynamic program over the same
recursion the enumerator uses, not by rejection over relation-length vectors --
and puts each through the cheap pipeline: named by the quipu theorem, carried by
the moves to one that is, or left over. The leftovers are the point, and each is
recorded with its overlap profile so their shapes can be counted afterwards.
Research H-019 is the hypothesis it is aimed at; the first 60 draws at `n = 12`
put the leftover rate at 28% against 0.6% at `n = 9`.

Every task writes an append-only ledger under `logs/`, one line per finished
draw, and resumes from it, so a run can be stopped and restarted and a run of
1000 can be widened to 2000 without redoing the first 1000. `--budget-hours`
stops cleanly and exits 2, which is what `overnight.py` restarts on.
`python batch.py --list` is the inventory of long jobs, including the ones that
keep their own front door.

### Searching without walking the same algebra twice

A mutation search enumerates mutation *sequences*, and many of them reach the
same algebra: at `n = 9` and depth 6 a walk visits 19483 nodes that are 1708
distinct algebras, and the ratio roughly doubles per level.

```python
from quivermutation import fingerprint, search

visited = fingerprint.Visited()
search.mutationSearchDepthFirst(algebra, 6, visited = visited, printOutput = False)
visited.summarise()          # nodes, distinct, skipped, ratio
```

The key is **exact**, not probabilistic, and the reason is that there is no
isomorphism problem: vertex labels do not move under mutation, so two algebras
reached from one start are equal on the nose or not at all. What is ambiguous is
the naming of parallel arrows -- a permutation within each bundle, and in
practice always one bundle of two -- and the **sign gauge**, rescaling an arrow
by `-1` being an automorphism that the procedure turns out to be sensitive to.
Both are quotiented out. It reaches exactly what the plain walk reaches, checked
over every LNA of `n = 5` to `8`, and is 6.1x faster at `n = 9`, depth 6.
`merges.py` uses it by default. See `research/` F-049, F-050 and E-042, E-043.

## Other families that could carry the classes the theorem misses

```bash
python families.py trees 9 10 11        # every tree, against every LNA
python families.py quipus 9             # every quipu with relations, against the
                                        # LNAs that lie in no quipu class
python families.py quipus 9 --verify 4  # and search each lead for a mutation path
python families.py members 9            # the quipu algebras a walk *proves* are in
                                        # each class the theorem does not name
python families.py free 8               # are two-arrow relations free on a quipu?
```

The quipu theorem names the class of an LNA with almost separate relations by a
**tree with no relations**. Everything it does not cover has to be classified
some other way, and `families.py` asks whether some other family of quivers
plays the same role for those. Both halves work by the Coxeter polynomial, which
is a derived invariant: a candidate whose polynomial no LNA of the length carries
is ruled out outright, and one that matches is a lead for a mutation search to
settle.

* `trees` -- every tree of the order, not only the quipus. Research F-031 did
  this for the trees of maximum degree three; this does the rest, and the answer
  is still no: **no tree outside the quipu shape carries the polynomial of an LNA
  outside a quipu class**, at orders 9 to 12 (F-033).
* `quipus` -- quipu quivers that **do** carry relations, in every orientation.
  Here the answer is yes, and in quantity: every Coxeter polynomial of an LNA in
  no quipu class is carried by quipu algebras with relations, thousands of them
  per polynomial (F-034).
* `members` walks the other way, out of the LNAs, so everything it prints is in
  the class it is printed under with a mutation path behind it. Every one of the
  9 LNAs outside a quipu class at `n = 9` and all 262 at `n = 10` reaches one
  within three mutations. Whether one of them is *canonical*, which is what would
  make this a theorem in the shape of the quipu theorem, is research H-014.

`free` is the caveat attached to that. `corollary:lengthtworelations` of
arXiv:2310.08346 makes a relation of two arrows free on a **line**; on a quipu
with a branch it is not, and about half the ideals that have one change their
Coxeter polynomial when it is deleted (F-035). So `--min-arrows 3`, which is what
makes order 11 affordable, is a real restriction of the family and not a
normalisation.

### Reorienting a tree is free

```python
from quivermutation import reflections as rf

rf.reflectionSequence(rf.arrowsOf(quiver), target)   # the mutations that reorient it
rf.mutationBridge(oneLNA, another, depth = 4)        # join two classes through their tree
```

A relation-free quiver on a tree mutates at a source by reversing exactly the
arrows there, which is the BGP reflection; left mutation at a sink is its
inverse. So the orientations of a tree are one **mutation** class, not only one
derived class -- which is what the classification needs, since it merges two
classes when both reach the same tree, and the two can reach it pointing
different ways. It does so in about a fifth of the pairs at n = 7 and n = 8, and
joining those orientations takes up to 11 mutations where a classification search
runs at depth 6. `mutationBridge` writes the whole path down. See research F-036.

### Where the Coxeter polynomial is not enough

```bash
python classify.py 9 --collisions
```

The Coxeter polynomial is a derived invariant but not a complete one, and it is
possible to say exactly where it fails without running any mutation. For a tree
it is determined by the tree's adjacency spectrum, so two cospectral
non-isomorphic quipus give algebras that are not derived equivalent yet share a
Coxeter polynomial. There are none below order 9 -- which is why the published
classification up to n = 8 can be read off the polynomial -- one pair at order 9,
two at order 10, four at order 11, thirteen at order 12.

The library underneath is usable directly:

```python
import quivermutation as qm
from quivermutation import nakayama as nk

a = nk.LinearNakayamaAlgebra(5, "300")      # 1->2->3->4->5, with 1->2->3->4 = 0
a.kupischSeries()                            # (3, 4, 3, 2, 1)
a.quipuName()                                # 'P^(2)_(1,1)', i.e. Dynkin D_5
a.coxeterPolynomial()                        # lambda**5 + lambda**4 + lambda + 1

qm.quiverMutationAtVertices(a, [4, 1])       # mutate, right at 4 then right at 1
table, report = qm.classifyLength(6)         # the whole classification
```
