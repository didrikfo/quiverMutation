# Findings

Established results, newest first. Each carries the evidence it rests on.
See [`README.md`](README.md) for conventions.

---

## F-015 — What the cheap steps cannot place is overlapping pairs, and no rule breaks one
*2026-09-16*

Seeding from the quipu theorem covers the LNAs whose consecutive relations
overlap in at most one arrow, so every row it misses carries an **overlapping
run**: a maximal group of relations chained by overlaps of two or more arrows.
Measuring the runs that actually occur says exactly what rule discovery should
be aimed at, and it turns out to be one shape.

    python unplaced.py 8
    python unplaced.py 9

| | n = 8 | n = 9 |
|---|---|---|
| LNAs | 429 | 1430 |
| placed by the theorem alone | 233 (54%) | 610 (43%) |
| placed by theorem + move orbits | 246 (57%) | 644 (45%) |
| still needing a search | 183 (43%) | 786 (55%) |

**Every unplaced row carries at least one overlapping run** — 182 of 183 carry
exactly one at n = 8, 772 of 786 at n = 9 — so the runs are not merely correlated
with the gap, they are the gap.

**One run dominates.** The commonest is the maximally overlapping pair of
length-3 relations, `(1:3) (2:3)`: it blocks 47 of the 183 unplaced rows at n = 8
(26%) and 159 of 786 at n = 9 (20%). The top 5 runs cover 54% of the unplaced
rows at n = 8 and 42% at n = 9; the top 10 cover 64% and 51%.

**The move table does reduce overlap, but almost never where it is needed.** Of
its 64 rules, 11 have a smaller maximum overlap after than before — so the
earlier guess that the moves simply preserve the obstruction is wrong. What is
true is narrower and sharper: **every one of those 11 needs either a third
relation in the window or relations of length 2**, and *none* of them applies to
a bare overlapping pair. Measured over the orbits: of the 820 LNAs of length 9
with an overlap of 2 or more, only **49 (6%)** reach a smaller overlap anywhere
in their orbit under the whole table.

That is why seeding plus move orbits adds only 3 points at n = 8 and 2 at n = 9
over seeding alone (H-003, now with numbers): the rules fire, they just do not
fire on the shape that is blocking the table.

**Reproduction.** `unplaced.py` for the table above; for the 11 rules and the
49, compare `discover.maximumOverlap` before and after each rule of
`lnaMoves.VERIFIED_MOVES`, and over `lnaMoves.orbitOf(9, ...)` for every LNA of
length 9.

---

## F-014 — n = 10: all 36 quipus found, and the predicted collision confirmed
*2026-09-14* · **partial — the merge step is unfinished**

`classifyLength(10)` placed all 4862 rows. After naming: **133 search classes**,
61 of them named — **all 36 quipus of order 10**, 13 classes of canonical type, 12
certified not piecewise hereditary. 72 classes remain unnamed and 17 groups
remain merge candidates, so this is **not yet a classification**; the resolve step
has not been run.

What is already established:

* every quipu of order 10 occurs as the class of some LNA of length 10, as H-005
  needs;
* the one group proved `separated` is `P^(1,4)_(1,0,2)` against
  `P^(3,3)_(1,0,1)` — **one of the two cospectral pairs F-010 predicts at order
  10**, arrived at by a completely different route;
* 13 canonical-type and 12 non-piecewise-hereditary classes, against 1 and 1 at
  n = 9, so the non-quipu part of the classification grows quickly (H-004).

The mismatch between the *two* predicted collision groups and the *two* reported
separations is what exposed R-007 — one of the reported ones was a false
separation.

---

## F-013 — The pair slide holds for every relation length, with two mutations
*2026-09-14*

Two relations of equal length `l` starting at consecutive vertices (maximally
overlapping) slide one arrow along the quiver, provided no other relation shares
an arrow with their span:

- two **right** mutations at the first relation's source move the pair one arrow
  **left**;
- two **left** mutations at the second relation's target move it one arrow
  **right**.

The mutation count is **two whatever `l` is**; only the window widens, to
`l + 2` arrows.

**Evidence.** `lnaMoves.verifyMove` at every window position of every LNA, for
`l = 2..7`, over four lengths each (`l+3` to `l+6`): 22 confirmations per
direction per length, no failures. Generated rather than listed, by
`lnaMoves.pairSlideRules`.

Note the direction: right mutations slide the pair **left**. See R-004.

Discovery found only `l = 3` and `l = 4` on its own, because those are the
window widths that fit at the lengths being searched — see H-001, which this
confirms.

---

## F-012 — Certificates of non-piecewise-heredity propagate by deleting vertices
*2026-09-13*

Corollary `removevertex` of [arXiv:2310.08346](literature/2310.08346-non-piecewise-hereditary-nakayama.md):
if a Nakayama algebra is piecewise hereditary, so is the algebra obtained by
deleting any one vertex. Contrapositively, if some one-vertex deletion is not
piecewise hereditary, neither is the algebra.

Reach, against the paper's two direct criteria alone:

| n | LNAs | direct | with deletion |
|---|---|---|---|
| ≤ 8 | | 0 | **0** |
| 9 | 1430 | 1 | 1 |
| 10 | 4862 | 22 | 24 |
| 11 | 16796 | 265 | **308** |

**Evidence.** `piecewiseHereditary.notPiecewiseHereditaryByDeletion` over every
LNA at each length. The zero below length 9 is the load-bearing check: the paper
states every LNA of length ≤ 8 is piecewise hereditary, so any certificate there
would mean the deletion construction drops or extends the wrong relation.

**What this does not give.** Deleting a vertex preserves piecewise heredity
*only* — not the derived equivalence class, the Coxeter polynomial, or the quipu.
The certificate says exactly "this class is not a quipu class", so it can
separate a class from every quipu class but can never merge two classes.

---

## F-011 — Every n = 9 class is named, and there are exactly 20
*2026-09-13*

The 1430 LNAs of length 9 fall into **20** derived equivalence classes:

- **18** quipu classes, one per quipu of order 9;
- **1** class (`3345000`, 8 members) of canonical type, tubular weight `(2,4,4)`;
- **1** class (`3033030`, 1 member) not piecewise hereditary at all.

**Evidence, three independent routes agreeing.**

1. `classifyLength(9)` assigns all 1430 rows and leaves nothing as a candidate.
2. The count: 18 quipus of order 9 exist; arXiv:2310.08346 states exactly one LNA
   of length 9 is not piecewise hereditary; the remaining class is of canonical
   type. 18 + 1 + 1 = 20.
3. **20 is exact, not a lower bound.** Two classes can only merge if they share a
   Coxeter polynomial. Across the 20 the polynomials are distinct except for the
   cospectral pair of F-010, which is provably two classes.

`3033030` is precisely the quiver `(**)` of arXiv:2310.08346 — relations
1→4, 3→6, 4→7, 6→9 — checked relation by relation.

---

## F-010 — The Coxeter polynomial fails exactly at cospectral quipus
*2026-09-13*

The Coxeter polynomial of the path algebra of a tree is determined by the tree's
adjacency spectrum. So two **cospectral** non-isomorphic quipus give algebras that
are not derived equivalent yet share a Coxeter polynomial — and among these
classes that is the *only* way it can fail.

| order | quipus | collision groups | quipus involved |
|---|---|---|---|
| ≤ 8 | 2..11 | **0** | 0 |
| 9 | 18 | 1 | 2 |
| 10 | 36 | 2 | 4 |
| 11 | 64 | 4 | 8 |
| 12 | 127 | 13 | 27 |
| 13 | 241 | 30 | 61 |

**Evidence.** `quipuForms.cospectralQuipuGroups(n)`. For orders 4 to 11 the
groups it finds are *exactly* the groups with equal Coxeter polynomials computed
the expensive way, through each algebra's Cartan matrix.

**This is why the published n ≤ 8 table is clean**: below order 9 there are no
cospectral quipus, so grouping by Coxeter polynomial is right there and nowhere
else. The smallest collision, at order 9:

    P^(1,4)_(1,0,1)  =  A_{9,(1,3)}^{(3,6)}   (class 3060000)
    P^(1,2)_(1,1,2)  =  A_{9,(1,4)}^{(3,4)}   (class 3004000)

Reproduce with `python classify.py 9 --collisions`.

---

## F-009 — The length-12 crash was cycles, not length
*2026-09-13*

`allRelsBetweenVertices` and `extendRel` recurse along the arrows out of a vertex
and tracked no visited set, so any cycle in the quiver made them descend forever.
`mutationSearchDepthFirst` then made that fatal rather than local by computing all
relations of a quiver at the top of every node *before* testing it for cycles, so
the first mutation producing a cyclic quiver killed the whole search on the next
node.

This explains the shape of the reported failure exactly: the crash appeared only
when running length 12, yet the offending quiver "didn't necessarily consist of
12 vertices". What mattered was that a mutation at that length finally produced
a cycle.

**Evidence.** `allRelsBetweenVertices` raises `RecursionError` on the three-cycle
`1→2→3→1` with `[1,2,3] = 0`. Fixed by bounding both recursions to simple paths
and testing for cycles before enumerating relations; the full n = 7 search is
byte-identical before and after, so nothing else moved. Pinned in
`tests/test_cycles.py`.

---

## F-008 — `reducePathAlgebra` preserves the algebra
*2026-09-13*

Reduction changes the quiver — it cancels arrows against inadmissible relations —
but must present the same algebra, and the Cartan matrix is the sharpest cheap
witness since it is indexed by vertices that do not move.

**Evidence.** Every legal mutation of depth ≤ 3 out of every LNA of lengths 5 to
8: **38095 reductions, zero changes** to the Cartan matrix, and the exact and
heuristic Cartan matrices agree on all of them.

---

## F-007 — Mutation reachability at bounded depth is directional
*2026-09-13*

`mutationSearchDepthFirst` walks only **right** mutations, so A can reach B at
depth `d` while B reaches nothing at that depth. Since
`rightMutate(dual(P)) = dual(leftMutate(P))` and the relation dual of an LNA is
derived equivalent to it, a right-mutation search out of `dual(X)` covers the
left-mutation directions out of `X`, and everything it reaches is in `X`'s class.

**Evidence.** At n = 8, seeding leaves class `340030` = A_{8,(1,2,5)}^{(3,4,3)}
reaching nothing classified at depth 7; searching from its relation dual finds
the link at depth 6 and closes the classification.

---

## F-006 — Two-path relations are differences, not sums
*2026-09-13*

A relation written as two paths means commutativity throughout this codebase —
`applyRelSetToPath` substitutes one for the other — so it means `p - q`, not
`p + q`.

Not cosmetic: three commutativity relations among three parallel paths read as
sums give `p = -q`, `r = -q` and `p + r = -2q`, forcing `q = 0` and collapsing a
Hom space that should be one-dimensional.

**Evidence.** Found via F-008: nine of 38095 reductions appeared to change the
Cartan matrix, and every one was this sign reading rather than a fault in the
reduction. See R-002.

---

## F-005 — Signs alone do not close; coefficients must be integers
*2026-09-13*

From `p + q + r = 0` and `p - q = 0` follows `2p + r = 0`, which has no expression
with coefficients in `{-1, 0, +1}`. The first step that combines two relations
leaves the sign-only world.

Integer coefficients are the smallest choice that closes and cost nothing over
signs — the same dict with a wider value type. The row reduction in
`relationAlgebra` is over the rationals already, so moving to a field would change
nothing.

---

## F-004 — The set-of-paths model misses zero relations through chained squares
*2026-09-13*

In the 2×2 commutative grid, the two commutativity relations make all three paths
from 1 to 6 equal, so adding `[1,2,3,6] = 0` kills all three. `pathHasZeroRel`
recognises only `[1,2,3,6]`, because it looks for a zero relation sitting
contiguously inside the path. `relationAlgebra.isInIdeal` gets all three.

A second symptom: `reducePathAlgebra` turns `{p,q,r}` plus `{p,q}` into `{p,q}`
and `{r}`, valid for `p+q+r=0, p+q=0` but not for `p+q+r=0, p-q=0` where the
answer is `2p+r=0`; `numberOfPathsUpToRels` reports 2 for that same algebra, so
the two halves of the code disagree about one object.

**Nothing published is affected**: see F-008.

---

## F-003 — The quipu theorem inverts, naming a class in O(1)
*2026-09-12*

Theorem `thm:QuipuToAn` of arXiv:2305.06642 read backwards: an LNA
`A_{n,(n_0..n_r)}^{(l_0..l_r)}` with almost separate relations comes from the
quipu with

    m_i = l_i - 2,  k_0 = n_0,  k_i = n_i - n_{i-1} - m_{i-1} - 1,
    k_{r+1} = n - n_r - m_r - 1

so it names its own class with no searching. Relations of length 2 are dropped
first, since they do not change the class.

**Evidence.** Round-trips on every quipu the paper names; constant on every
published class; agrees with what a mutation search reaches wherever the search
can reach an answer. Coverage falls with length (100% of rows at n = 4, 54% at
n = 8, 19% at n = 12) but the set of quipus it names does not: it finds every
class at every length checked.

---

## F-002 — Step 3's cyclic case is not implemented
*2026-09-12*

For a minimal relation `r: i ⇢ i` the procedure calls for one arrow
`(α r̄): i* → t(α)` per arrow `α` out of `i`; the code adds a single arrow from
`r`'s source to its target, which on a cycle is a loop on `i*`. The second worked
example of arXiv:2112.08129 does not reproduce.

Pinned as a strict xfail in `tests/test_mutation_procedure.py`. The LNA search
never reaches it, since it stops descending at the first cycle, so no published
result depends on it.

---

## F-001 — A NetworkX change had silently broken every Coxeter polynomial
*2026-09-12*

Since NetworkX 3.1, `all_simple_paths(G, v, v)` yields the trivial path `[v]`
where it used to yield nothing. `numberOfPathsUpToRels` counted it, putting 2
instead of 1 on the Cartan diagonal, so `det(C) = 2^n` instead of 1 and every
Coxeter polynomial came out with fractional coefficients.

**Evidence.** After the fix, A_5 gives `1 + λ + … + λ^5` and D_5 gives
`(λ+1)(λ^4+1)`, as they must.
