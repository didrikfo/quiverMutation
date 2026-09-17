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
