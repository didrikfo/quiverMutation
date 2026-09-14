# Findings

Established results, newest first. Each carries the evidence it rests on.
See [`README.md`](README.md) for conventions.

---

## F-015 — The procedure on coefficients, and the two relations the old one missed
*2026-09-14*

`procedure.py` is steps 1-7 of arXiv:2112.08129 and the cleanup after step 7,
with relations as `relationAlgebra` combinations instead of sets of paths. The
coefficients come out of the steps rather than being guessed back, and step 7 is
computed as the kernel the paper says it is. It is what
`mutation.quiverMutationAtVertex` and `reduction.reducePathAlgebra` now run.

**Agreement with the implementation it replaced**, gated on
`mutationIsPossibleAtVertex` throughout so both walk the same mutations:

| check | cases | differences |
|---|---|---|
| one mutation at every admissible vertex of every LNA, n = 4..8 | 2349 | 0 |
| depth-3 walks, n = 5 | 1446 steps | 0 |
| depth-3 walks, n = 6 | 7496 steps | 0 |
| depth-3 walks, n = 7 | 37470 steps | **2** |
| the exact cleanup on the old steps' output, n = 5,6 depth 3 and n = 7 depth 2 | 13652 | 0 |

**The two differences are relations the old implementation missed**, and the
paper's step 7 says so — see R-007. Both are at n = 7 after three mutations:

    A_{7,(1,4)}^{(4,3)} = 40030, mutated at 1, 4, 2
      old: 1;2;5;6;7 | 3;2;5;4 | 5;4;7 = 5;6;7
      new:   2;5;6;7 | 3;2;5;4 | 5;4;7 = 5;6;7     (and 2;5;6;7 implies the old one)

    A_{7,(1,4,4)}^{(4,3)}-ish = 44030, mutated at 1, 4, 2
      new has 2;5;6;7 = 0 in addition to everything the old one has

In the first, `2;5;6;7` is **not** in the ideal the old answer generates, so the
two are different algebras: total dimension 23 against 22, where the algebra
being mutated has 22.

**Neither is caught by the checks that were in place.** Both preserve the
Coxeter polynomial, both survive a there-and-back left mutation, and both reach
`P^(1,2)_(1,0,1)` — the quipu the theorem names for that class — by a depth-7
search. Which is worth recording on its own: *the Coxeter polynomial, the round
trip and the hereditary form can all three be satisfied by an algebra that is
not the mutation*, because all three see only the derived equivalence class, and
a missing relation can leave the class unchanged.

**The coefficients never change a Cartan matrix.** Computing it from the true
combinations and from `fromPathSet`'s guess over the same algebra agrees on every
quiver reached within depth 3 of every LNA of lengths 5 and 6 — 1239 of them,
zero differences. So no published Coxeter polynomial was ever wrong because of
the guess; what the guess cost was the reasoning, as R-003 records.

**It is faster, not slower.** At n = 7 over the 462 admissible single mutations:
the procedure 0.26 s against 0.94 s, the admissibility condition 0.14 s against
1.74 s. Exact linear algebra over the ideal beats the hand-rolled list surgery it
replaces, by 3.6x and 12x. That was the opposite of what was expected, and is
why the switch was affordable.

**The classification is unchanged.** `python classify.py 8` on the new procedure
gives the same 11 classes with the same sizes --
133 + 65 + 64 + 64 + 40 + 26 + 13 + 10 + 9 + 4 + 1 = 429 -- as the partition
recorded in NOTES.md for the old one, and n = 7 gives the same
54 + 32 + 29 + 7 + 6 + 4 = 132. So the two missing relations did not change a
published class: they were lost on quivers the search passes through, not on the
LNAs it records. n = 6 now takes 5 seconds against 13, and n = 7 40 seconds
against 37.

**Evidence.** `tests/test_procedure.py` (agreement at n = 4..6, with n = 7 and 8
marked slow; step 5's minus sign; step 7 as a kernel; the cleanup; the two
recovered relations). The wider runs are E-014.

---

## F-014 — The quipu's end exchanges are already exactly right, and they do not merge the n = 9 pair
*2026-09-14*

The quipu notation is not unique, and one of the re-readings is the one that
matters most here: at the **outermost** foot the main-string end segment and the
cord are the only two branches, so exchanging them — `k_0` with `m_0`, or
`k_{r+1}` with `m_r` — is an isomorphism of the tree and the two LNAs it names
are derived equivalent. **This is not a symmetry the code was missing.** It was
already subsumed by canonicalisation, and it is now checked directly, from both
sides:

1. **Canonicalisation is exactly tree isomorphism.** Over *every* parameter pair
   of each order 3–11 — 4180 names at order 9, 28656 at order 11 —
   `quipuForms.quipuParameters` gives two names the same parameters **if and only
   if** `networkx` says their graphs are isomorphic. Both directions, no
   exceptions. So no end exchange (nor the backwards reading, nor any other
   relabelling) can split a class, and no pair of genuinely different trees is
   being run together.
2. **The LNA-side operations agree with it exactly.** `cor:EquivNakayamaAlgebras`
   is now implemented directly on the relations, independently of `quipuForms`,
   as `LinearNakayamaAlgebra.swapFirstRelation`, `swapLastRelation`,
   `relationDual` and `withoutShortRelations`, with `classPreservingOrbit` for
   their closure. For every LNA of length 4–10 with almost separate relations and
   no length-2 relation, **the orbit of the paper's operations equals the set of
   LNAs that `quipu` names with the same quipu** — 64 algebras in 18 orbits at
   n = 9, 128 in 36 at n = 10, and the largest orbit is 8, which is the paper's
   own bound.

The exchange is at the ends only. Applied to an interior gap `k_i`,
`0 < i < r + 1`, it changes the tree, because the foot there has a third branch
running on along the main string: `P^(1,2,1)_(1,1,2,1)` exchanged at `k_1` is a
different tree of the same order.

**Consequence for the n = 9 pair.** The two classes R-006 says the workbook
over-merged are not related by any of this:

    P^(1,4)_(1,0,1)   (3060000)   diameter 6   -- exchanges to P^(1,1)_(1,0,4)
    P^(1,2)_(1,1,2)   (3004000)   diameter 5   -- both exchanges fix it

The end exchange at the last foot of the first one does fire, and lands on
`P^(1,1)_(1,0,4)` = `3030000`, which canonicalises straight back to
`P^(1,4)_(1,0,1)` — so it adds `3030000` and `0003030` to that class, where they
already are. It does not reach the other quipu, and nothing can: the trees have
different diameters, so they are not isomorphic, and for hereditary algebras of
tree type the underlying tree *is* the derived equivalence class.

The same split is **already in the published table** one vertex down:
`P_(1,0,3)^(1,1)` and `P_(1,1,2)^(1,1)` are two separate rows of the n <= 8
classification in arXiv:2305.06642, and adding a vertex to the last cord of each
gives exactly this pair. n = 9 is not a new claim — it is the published one at
the first order where the Coxeter polynomial can no longer see it (F-010).

**Evidence.** `tests/test_quipu_symmetry.py`, 44 tests at orders/lengths up to 9
plus three marked slow at 10 and 11. Reproduce with

    .venv/bin/python -m pytest tests/test_quipu_symmetry.py -q

**What this does not settle.** Both routes above descend from `thm:QuipuToAn`
being correctly inverted, and the third route — reaching a relation-free quiver
by mutation — is out of range for this pair: searches from all four long-relation
members of `P^(1,4)_(1,0,1)` and both of `P^(1,2)_(1,1,2)` reach no hereditary
quiver at depth 6 (E-013). A genuinely theorem-free separation would need a
derived invariant computed from the algebra itself; the Avella-Alaminos–Geiss
invariant for gentle algebras is the candidate (NOTES idea 22).

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

*Amended 2026-09-14.* The evidence above stands, but the **regression test** for
it did not: `test_each_rule_holds_wherever_it_applies` asked for lengths 5 to 8
whatever the rule, which is `l+3 .. l+6` only for `l = 2`. The window of the
pair slide is `l + 2` arrows wide and needs `l + 3` vertices to sit in, so the
rules for `l = 6..9` got **zero confirmations** at those lengths and the test
failed on all eight of them — for asking where they cannot occur, not for
anything wrong with the rules. It had been failing since the family was
generated up to `l = 9`; it is in the slow set, which is why it went unseen. The
test now derives its lengths from each rule's own width, and a second test
asserts that it always can.

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
