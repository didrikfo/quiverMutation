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

## Running a mutation class search

`mutationSearch` walks every LNA of a given length, runs a depth-first search of
tilting mutations from it, and records every other LNA reached in a CSV file.
It writes its output to the current working directory, so run it from a scratch
directory:

```bash
mkdir -p out && cd out
python -c "import quiverMutation as qm; qm.mutationSearch(6, 6, 0, createNewCSVfile=True)"
```

This produces `A_6_mutation_classes.csv` with one row per LNA, giving the
mutation class it was assigned to, the mutation path from the class
representative, its Coxeter polynomial, and the vertex numbering.

Classes sharing a Coxeter polynomial are candidates for merging: the Coxeter
polynomial is a derived invariant but not a complete one, and the depth-first
search is not guaranteed to find every mutation path between two LNAs of the
same class. See `NOTES.md` for the current state of that workflow.
