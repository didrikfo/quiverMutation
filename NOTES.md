# Working notes

State of the codebase, what the model can and cannot express, and the backlog.
Written while picking the project back up; keep it current.

## What the code does

The repo implements the combinatorial rule for tilting mutation of
[arXiv:2112.08129](https://arxiv.org/abs/2112.08129), and applies it to classify
linearly oriented Nakayama algebras (LNAs) up to derived equivalence, which is
how the data behind [arXiv:2305.06642](https://arxiv.org/abs/2305.06642) was
produced.

### The model of a path algebra

`pathAlgebraClass.PathAlgebra` holds

* `quiver`: a `networkx.MultiDiGraph` whose nodes are the vertices (integers,
  `1..n` for a line);
* `rels`: a list of relations.

A **relation** is a list of paths, and a **path** is a list of vertices. So
`[[1,2,4],[1,3,4]]` is the commutativity relation between the two paths from 1
to 4, and `[[1,2,3]]` is the zero relation on the single path `1 -> 2 -> 3`.
By convention `rel[0][0]` is the source and `rel[0][-1]` the target, and every
path in a relation shares both.

This is the deliberate simplification: a relation is a **set** of paths, not a
linear combination. There are no coefficients, so the model cannot tell
`p - q = 0` from `p + q = 0`, and cannot express `2p - 3q + r = 0` at all. A
relation with one path means that path is zero; a relation with two or more
means they are identified up to sign. See "Improving the model" below.

A vertex has an implied identity path, which nothing in `rels` represents; the
`+1` on the diagonal of `cartanMatrix` stands in for it.

### The mutation procedure

`quiverMutationAtVertex(pathAlg, vertex)` applies steps 1-7 of the paper's
procedure and returns a new `PathAlgebra`. It does not clean up after itself:
the result routinely contains inadmissible relations (ones containing a path of
length one) and redundant relations.

`reducePathAlgebra(pathAlg)` does the cleanup, in this order: drop illegal
relations, dedupe relation paths, cancel each length-one path in a relation
against the corresponding arrow (the "Note" after the paper's step 7),
substitute through the arrows that cancelled, drop non-minimal zero relations,
drop redundant relations and existing subrelations, iterate to a fixed point.

`quiverMutationAtVertices(pathAlg, vertices)` mutates at each vertex in turn,
reducing after each. A negative entry means left mutation, which
`leftQuiverMutationAtVertex` performs as dual -> right mutation -> dual.

`mutationIsPossibleAtVertex(pathAlg, vertex)` is the admissibility test: an
arrow out of the vertex must exist, the quiver must have no parallel arrows, and
`Hom(P_i*[1], Lambda)` must vanish.

### The algebra classes

`nakayama.LinearNakayamaAlgebra(length, relLengths)` is a `PathAlgebra` that
knows it is an LNA. It owns the three names that were being converted between by
hand all over the module -- the per-vertex relation lengths `[2,2,3,0,0]`, the
class name `'22300'`, and the relation string `'1;2;3|2;3;4|3;4;5;6'` -- along
with the Kupisch series, the relation dual, whether the relations are almost
separate, the quipu, and the Cartan matrix and Coxeter polynomial. It is
hashable and compares by structure, so LNAs work as dict keys.

`nakayama.QuipuAlgebra(k, m)` is the path algebra of the quipu quiver
`P^(m)_(k)`, with no relations. `QuipuAlgebra.fromLNA` and `correspondingLNA`
are the two directions of the theorem.

Checked: for every LNA of length <= 7 with almost separate relations, the LNA
and its quipu algebra have the same Coxeter polynomial, computed from their own
Cartan matrices by entirely separate routes, and the quipu round-trips.

### The LNA search

`mutationSearch(lineLength, mutationDepthStart, startRow, createNewCSVfile)` is
the entry point, and it returns a `MutationClassTable`.

1. `generateAllPossibleLineRelations(n)` enumerates every admissible ideal on the
   linear quiver `1 -> ... -> n`. There are Catalan(n-1) of them.
2. `createMutationClassCSV(n)` writes one row per LNA to
   `A_<n>_mutation_classes.csv`, with empty class columns.
3. For each row still without a class, `mutationSearchDepthFirst` walks
   mutations to the given depth and collects every LNA it reaches into the list
   passed as its `collected` argument.
4. `mutationListLineCleanup` normalises each of those back to the standard
   numbering `1..n` and dedupes.
5. `assignMutationClassInTable` writes the class name, the mutation path from
   the class representative, the Coxeter polynomial and the vertex numbering
   into every row it reached.

`mutationClassTable.MutationClassTable` owns the table: one row per LNA keyed by
its relation string, with a dict index over that key and a polars DataFrame for
interchange (`toDataFrame`, `writeParquet`). `classesByCoxeterPolynomial` is the
entry point for the merge step.

A class is named after the LNA that seeded it, in the per-vertex relation-length
notation: `A7_22300` is the line on 7 vertices with relations of 2, 2 and 3
arrows starting at vertices 1, 2 and 3. The same notation appears as the
`relationList` argument of `lineQuiverExample`.

The depth-first search is bounded by `depth`, which `mutationSearch` decays as
`max(depthStart - floor(log10(row)), 2)` to keep later rows affordable, and it
stops descending as soon as a cycle appears in the quiver.

Because the depth is bounded, **the search finds a lower bound on each class**:
two LNAs in the same class can end up in different search classes if no mutation
path between them fits in the depth. Merging those was the manual step. The
Coxeter polynomial identifies the candidates, since it is invariant under
derived equivalence — but it is not a complete invariant, so agreement is not
proof. It does happen to separate every class for n <= 8.

## Verified against the papers

* `generateAllPossibleLineRelations(n)` returns Catalan(n-1) relation sets for
  n = 2..11.
* The acyclic worked example of arXiv:2112.08129 (mutate the 7-vertex quiver at
  vertex 3) reproduces the paper's quiver and relations exactly.
* The Coxeter polynomial does not move along any legal mutation path of depth
  <= 3 out of 14 different LNAs of length 5 to 7.
* Every class in the published n <= 8 table has a single Coxeter polynomial, and
  distinct classes have distinct ones.
* n = 5: the search finds both classes outright, 8 + 6 = 14 LNAs.
* n = 6: the search finds 5 classes with 4 distinct Coxeter polynomials. Merging
  the two that share the D_6 polynomial gives the paper's 4 classes for n = 6
  (A_6, D_6, E_6, D~_5) with sizes 16 + 13 + 12 + 1 = 42.
* n = 7: the search finds 11 classes with 6 distinct Coxeter polynomials,
  matching the 6 quipus of order 7 in the paper, sizes
  32 + 29 + 54 + 7 + 6 + 4 = 132.
* n = 8: the search finds 28 classes with 11 distinct Coxeter polynomials,
  matching the 11 quipus of order 8 in the paper. Every one of the paper's
  classes lands inside a single group, and the group sizes
  1 + 4 + 9 + 10 + 13 + 26 + 40 + 64 + 64 + 65 + 133 account for all 429 LNAs.
  The whole run takes about 7 minutes at depth 6.

Run times at depth 6, after the deepcopy fix: n = 5 a few seconds, n = 6 about
20 seconds, n = 7 90 seconds, n = 8 7 minutes. Each is pinned as a `slow` test.

## Known gaps and limitations

* **Step 3's cyclic case is not implemented.** For a minimal relation
  `r: i --> i` the procedure calls for one arrow `(alpha r-bar): i* -> t(alpha)`
  per arrow `alpha` out of `i`; the code adds a single arrow from `r`'s source to
  its target, which on a cycle is a loop on `i*`. Pinned as a strict xfail in
  `tests/test_mutation_procedure.py`. The LNA search never reaches it because it
  stops descending at the first cycle, so this has never affected a published
  result.
* **`mutationIsPossibleAtVertex` is stricter than the paper** when the vertex has
  more than one arrow out of it: it rejects the vertex as soon as one arrow out
  of it kills a nonzero path, where the paper only asks that some arrow out of it
  does not. The two agree when the vertex has a single arrow out of it, which is
  the only case the LNA search meets.
* **Parallel arrows are rejected outright**, because a path is a vertex sequence
  and so cannot name which of two parallel arrows it uses.
* **The recursive loop at length 12.** A full n = 12 run crashed in a recursive
  loop, somewhere in the mutation of relations. The offending LNA is not
  recorded. Reproducing it is the first step; a long-running n = 12 search with
  the recursion limit lowered and the failing quiver dumped on `RecursionError`
  would find it.
* **Relations are compared by sorted vertex sequences**, so two relations that
  are equal as ideals but written differently are distinct objects. Several
  functions exist mainly to paper over this (`removeDuplicateRels`,
  `removeDuplicateRelPaths`, `removeRedundantRelations`,
  `removeExistingSubrelations`, `minimizeCommutingRelation`).
* **The older entry points still round-trip through text files**, parsed by
  string slicing at fixed offsets in `readMutationsFromFile`.
  `mutationSearch` no longer does -- it collects in memory -- but
  `findMutationClassesForLine`, `collectMutationClasses`,
  `combineLineMutationFiles` and the scratch code in `main.py` still do.
* **`combineMutationClassesInCSVfile` does not work.** It was the start of an
  automated merge step and was never finished: it has a `#wrong!` marked append,
  a call to `quiverMutationAtVertices` missing its second argument, and a
  `baseClass.split('')` that raises. Treat it as a sketch of idea 12, not as
  code to fix.

## Backlog

Roughly in the order that unblocks the most.

### Performance and correctness

1. **Memoize the path enumeration.** With the graph deepcopy gone, the hot spot
   is `nx.all_simple_paths`, called 93k times in one depth-4 search on quivers
   that have not changed between calls. Needs a cheap structural key for a path
   algebra — which is also what a proper `__hash__`/`__eq__` on `PathAlgebra`
   would give.
2. **Find and fix the length-12 recursive loop**, per above.
3. **Canonical form for a path algebra.** A normalised, hashable representation
   would replace the dedupe-by-list-comparison machinery, let the search
   memoize on quivers it has already visited, and make `reducePathAlgebra`
   testable as "reduces to the canonical form".
4. **Avoid re-deriving relations from scratch after each mutation.** The search
   recomputes every relation between every pair of vertices at each node; most
   of that is unchanged by a mutation at a single vertex.

### The model

5. **Linear combinations in relations.** Half done. `relationAlgebra` has the
   value type and the exact ideal arithmetic (see "Coefficients" below); what
   remains is rewiring the mutation procedure and `reducePathAlgebra` to use it
   instead of the set-of-paths model.
6. **Name arrows.** Paths as vertex sequences cannot express parallel arrows or
   distinguish two arrows with the same endpoints. Arrow identities would lift
   that restriction and make the quiver a plain `DiGraph` of named arrows.

### OOP and structure

7. **Make `PathAlgebra` carry its own behaviour.** It is currently a thin
   container and every operation is a free function taking it as the first
   argument. Mutation, reduction, invariants and admissibility are all methods.
8. **Subclass for LNAs**, holding `n` and the per-vertex relation lengths, with
   the standard numbering, the relation-string form, the Kupisch series and the
   relation dual as its own operations. Most of the `...Line...` free functions
   collapse into it.
9. **Subclass for quipu quivers**, holding the `P^{(m)}_{(k)}` parameters, the
   orientation, and the CR-swap of arXiv:2305.06642 section "Cord/relation-swap"
   as a method. This is what makes idea 12 possible.
10. **Split the 2600-line module** along the seams that already exist: the
    procedure, relation algebra, invariants, the line search, quipus, IO.
11. **Replace the text-file round trip** in the remaining entry points
    (`findMutationClassesForLine`, `collectMutationClasses`,
    `combineLineMutationFiles`) the way `mutationSearch` now does it, with a
    `collected` list rather than a transcript parsed back by string slicing.

### Coefficients

`relationAlgebra` models a relation as `dict[path, int]` -- a linear combination
of paths with integer coefficients -- and decides ideal membership exactly, by
linear algebra rather than by pattern matching.

**Signs alone do not close.** From `p + q + r = 0` and `p - q = 0` follows
`2p + r = 0`, which cannot be written with coefficients in `{-1, 0, +1}`. The
first step that combines two relations leaves the sign-only world, so integer
coefficients are the smallest choice that works; they cost nothing over signs,
being the same dict with a wider value type. Rationals would serve equally and
the row reduction is over the rationals already.

**What the set-of-paths model gets wrong.** In the 2x2 commutative grid

    1 -> 2 -> 3          relations   [1,2,5] = [1,4,5]
    |    |    |                      [2,3,6] = [2,5,6]
    v    v    v
    4 -> 5 -> 6

all three paths from 1 to 6 are equal, so adding `[1,2,3,6] = 0` kills all
three. `pathHasZeroRel` recognises only `[1,2,3,6]`, since it looks for a zero
relation sitting contiguously inside the path.
`relationAlgebra.isInIdeal` gets all three. This is the failure described as
"a long zero relation which passes through multiple commutativity relations".

A second symptom, with three parallel paths and a mix of relation orders: the
current `reducePathAlgebra` turns `{p,q,r}` together with `{p,q}` into `{p,q}`
and `{r}`, which is valid for `p+q+r=0, p+q=0` but not for `p+q+r=0, p-q=0`,
where it should give `2p+r=0`. `numberOfPathsUpToRels` meanwhile reports 2 for
that algebra, so the two halves of the code disagree about the same object.

**Reading a relation without coefficients.** `fromPathSet` has to guess the
signs, and the guess matters:

* one path -- that path is zero, no sign to choose;
* two paths -- a **difference**, `p - q = 0`. This is what the rest of the repo
  means: `applyRelSetToPath` substitutes one path for the other, and
  `numberOfPathsUpToRels` treats a two-path relation as identifying them.
  Reading it as a sum is not harmless. Three commutativity relations among three
  parallel paths `p, q, r` become `p = -q`, `r = -q` and `p + r = -2q`, which
  forces `q = 0` and collapses a Hom space that should be one-dimensional. This
  showed up as nine apparent failures in the check below, every one of them the
  sign reading rather than a fault in the code being checked;
* three or more paths -- a **sum**, which is what step 4 produces. The true signs
  are not recoverable, and that is the central reason to move the procedure onto
  coefficients.

**How much of this matters for the published results: none of it so far.**
`relationAlgebra.cartanMatrixExact` agrees with the existing `cartanMatrix` on
all 624 LNAs of length <= 8, and on all 8101 quivers reached by walking every
legal mutation path of depth <= 3 out of all 188 LNAs of length 5 to 7. So no
Coxeter polynomial in the tables moves. The shapes where the two differ have not
turned up in an LNA search yet -- consistent with the crash only appearing at
length 12.

### Classifying a length

`classifyLength(n)` runs the whole thing and returns `(table, report)`:

1. **Seed.** `seedTableFromQuipuTheorem` fills in every LNA with almost separate
   relations, straight from the theorem, before a single mutation is computed.
   Each class is named by its quipu rather than by an arbitrary LNA, so two
   seeded classes with the same quipu are literally the same class. Coverage
   falls with length -- 100% at n=4, 54% at n=8, 19% at n=12 -- but the set of
   quipus it names does not: it already finds every class at every length
   checked.
2. **Search.** `mutationSearch` handles only the rows left over. Each inherits a
   seeded class the moment its search reaches a seeded LNA.
3. **Name.** `annotateHereditaryForms` gives a form to any class the search
   created that the theorem did not cover.
4. **Resolve.** `resolveMergeCandidates` takes each group of classes sharing a
   Coxeter polynomial that the hereditary form does not settle, and searches
   deeper for a mutation path between them.

**Reachability is directional**, and this matters. `mutationSearchDepthFirst`
walks only *right* mutations, so A can reach B at depth d while B reaches nothing
at that depth. Since `rightMutate(dual(P)) = dual(leftMutate(P))` and the
relation dual of an LNA is derived equivalent to it, searching from the dual
covers the other direction; `resolveMergeCandidates` does both. At n = 8 this is
exactly what settles the last class: `340030` = A_{8,(1,2,5)}^{(3,4,3)}, whose
relations overlap too much for the theorem, reaches nothing seeded, but its dual
finds the link at depth 6.

Results, all matching the published table with nothing left as a candidate:

| n | LNAs | classes | class sizes | time |
|---|------|---------|-------------|------|
| 6 | 42   | 4       | 16, 13, 12, 1 | ~13s |
| 7 | 132  | 6       | 54, 32, 29, 7, 6, 4 | ~37s |
| 8 | 429  | 11      | 133, 65, 64, 64, 40, 26, 13, 10, 9, 4, 1 | ~4min |

### The hereditary form

When a search leaves a quiver with **no relations**, the algebra is hereditary,
and for tree-shaped quivers derived equivalence is settled: two path algebras of
trees are derived equivalent exactly when the trees are isomorphic as undirected
graphs, since their orientations are related by BGP reflections. So the
underlying graph of any relation-free quiver a search reaches is a **complete**
derived invariant of the class, where the Coxeter polynomial is only a necessary
condition.

`quipuForms.canonicalTreeForm` encodes a tree canonically (AHU, rooted at the
centre, smaller of the two encodings when there are two centres), and
`quipuForms.quipuParameters` recovers the paper's `P^(m)_(k)` notation where the
tree is a quipu, canonicalised by taking the lexicographically smallest
parameter pair over every valid reading of the main string, since the notation
does not determine the quipu. `graphFromQuipuParameters` goes back, so a quipu
named in the paper can be compared with one a search found.

`mutationSearchDepthFirst` collects these for free during a class search via its
`collectedHereditary` argument, and `annotateHereditaryForms` fills in the
classes that search missed by iterative deepening from each member in turn.
`mergeReport` then turns the table into a decision:

* `certain` — classes reaching the same hereditary form. Provably one class:
  the search just missed the mutation path. **Merge these.**
* `separated` — classes sharing a Coxeter polynomial but reaching different
  hereditary forms. Provably distinct. **Do not merge these**, whatever the
  polynomial says. This is the case that makes the polynomial an incomplete
  invariant.
* `candidate` — classes sharing a polynomial where at least one has no
  hereditary form yet. Still needs work.

For n = 6 this reduces the entire hand-merge step to one `certain` entry and
nothing else, agreeing with the published table.

### Features

12. ~~**Better CSV post-processing, to cut the manual work.**~~ Done, by the
    hereditary form rather than by a deeper search: see "The hereditary form"
    below. `mergeReport` splits the same-polynomial groups into `certain`
    (provably one class), `separated` (provably distinct) and `candidate` (not
    yet settled). What remains is to shrink `candidate` — a class with no
    hereditary form reached needs either a deeper targeted search or the quipu
    seeding of idea 13.
13. ~~**Seed the table from quipu quivers.**~~ Done: `seedTableFromQuipuTheorem`,
    used by `classifyLength`. See "Classifying a length" below.
14. **More invariants, to separate classes the Coxeter polynomial cannot.**
    Candidates: the determinant and elementary divisors of the Cartan matrix,
    the Euler form, the number of indecomposables / the shape of the AR quiver,
    Hochschild cohomology dimensions, the derived invariants of Avella-Alaminos
    and Geiss for gentle algebras (LNAs are gentle), and the silting/tilting
    quiver's local structure. Each wants to be a function
    `PathAlgebra -> hashable`, so the merge step can key on a tuple of them.
15. **Certificates both ways.** A merge should record the mutation path that
    proves the equivalence; a split should record the invariant that separates
    the two classes. Then a classification is checkable without rerunning it.
16. **Push past n = 11.** Catalan growth means n = 12 is 58786 LNAs and n = 15
    is 2674440, so the search has to get cheaper per LNA and the table has to
    stop being a CSV read into memory.
