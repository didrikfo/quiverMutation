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

The classification runs in four steps, described in `NOTES.md`: seed every LNA
the quipu theorem of arXiv:2305.06642 covers, search by mutation for the rest,
name any class the theorem missed by the hereditary algebra its search reaches,
and settle whatever is left by a deeper search. For n <= 8 this reproduces the
published classification with nothing left over, replacing what used to be a
hand-merge over the CSV.

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
import nakayama as nk
import quiverMutation as qm

a = nk.LinearNakayamaAlgebra(5, "300")      # 1->2->3->4->5, with 1->2->3->4 = 0
a.kupischSeries()                            # (3, 4, 3, 2, 1)
a.quipuName()                                # 'P^(2)_(1,1)', i.e. Dynkin D_5
a.coxeterPolynomial()                        # lambda**5 + lambda**4 + lambda + 1

qm.quiverMutationAtVertices(a, [4, 1])       # mutate, right at 4 then right at 1
table, report = qm.classifyLength(6)         # the whole classification
```
