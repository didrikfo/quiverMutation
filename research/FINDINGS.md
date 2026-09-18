# Findings

Established results, newest first. Each carries the evidence it rests on.
See [`README.md`](README.md) for conventions.

---

## F-046 — The literature merges two pairs at `n = 11` that no move we have reaches
*2026-09-19*

The point of the sweep, and the only part of it that hands us something the
pipeline does not already do. Three published equivalences between radical-power
LNAs at `n = 11`:

| | source |
|---|---|
| `N_11(3) ≃ N_11(7)` | arXiv:2112.15587 Prop. 4.1(2), `a = 3`, `b = 7`, `n = (a-1)(b-1) = 12`, taken at `n - 1` |
| `N_11(6) ≃ N_11(7)` | arXiv:2203.15735 Prop. 4.5, `N(2r-1, r) ≃ N(2r-1, r+1)` at `r = 6` |
| `N_11(4) ≃ N_11(5)` | arXiv:2112.15587 Prop. 4.1(2), `a = 4`, `b = 5`, `n = 12`, at `n - 1` |

In relation-length lists (**not** digit strings — a relation of ten arrows has no
digit, and this is `n = 11`):

    N_11(3) = (3,3,3,3,3,3,3,3,0)      N_11(6) = (6,6,6,6,6,0,0,0,0)
    N_11(4) = (4,4,4,4,4,4,4,0,0)      N_11(7) = (7,7,7,7,0,0,0,0,0)
    N_11(5) = (5,5,5,5,5,0,0,0,0)

**Against the full union-find partition of `n = 11`** — the rule table, the edge
moves, the double mutation and the free move, which is every move we have —
`N_11(3)` and `N_11(6)` share an orbit of 954 rows, and `N_11(7)` sits in a
different orbit of 406, `N_11(4)` in one of 461 and `N_11(5)` in one of 74. So
the three equivalences above **merge four of our orbits into two**, and not one of
those merges is a move.

This is not a budget artefact. `derivedOrbits` is a union-find over every LNA of
the length, not a bounded walk, and it is the same computation the classification
runs on. As a second reading, `movesJoin` at a cap of 400000 rows **exhausts**
both orbits for `N_11(3)` against `N_11(7)` — 948 and 34 rows walked forward,
cap never approached — so the walk ran out of moves, not of budget, which is the
distinction E-037 section 3 was burned by.

**Nothing else we have reaches them either.**

* Coxeter polynomials agree in all three pairs, as they must, so the invariant is
  silent rather than confirming.
* **No tree of order 11 carries that polynomial**, so these classes are outside
  the quipu theorem entirely — this is the region the project is actually stuck
  in, not a corner of it.
* All five algebras are derived **wild** (Euler form indefinite), so F-045's
  criterion says nothing about them by construction.
* The quipu theorem does not name them, `canonicalWeightType` does not, and the
  moves do not.

**So: at `n = 11`, three lines of a paper do what the whole pipeline cannot.**
That is worth stating plainly, because every other result of this sweep runs the
other way — F-043 found the literature's radical-power merge already inside the
move table, and F-044 found its classification already reproduced by our naming.
Here the direction reverses.

**What to do with it.** Two things, in order.

1. *Use them.* Three merges is three merges; fold them into the seeding the way
   F-032's double mutation was folded in. They are cheap: both sides are named by
   a closed form in `(n, r)`, so the rule is a lookup, not a search.
2. *Ask what move they would be.* Each pair is two rows of the same length whose
   relation lengths differ throughout, joined by a proof that goes through
   `vect-X(2,a,b)` and one-point extensions rather than through mutation. If
   there is a move behind `N_n(a) ≃ N_n(b)` for `n = (a-1)(b-1) ± 1`, the search
   has never found it, and a targeted search between exactly these pairs — which
   is what `meetingPoints` is for — is the cheapest way to ask.

Reproduce: `hs.py` and `hsgap.py` in the sweep's scratchpad; the authoritative
line is `freeMoves.derivedOrbits(11, free = True, edges = True, doubles = True)`
and comparing the two orbit keys. E-040,
`research/literature/2112.15587-nakayama-fuchsian-singularities.md`,
`research/literature/2203.15735-one-branch-extensions-rectangles.md`.

---

## F-045 — Brüstle's invariant classifies the derived-tame LNAs outright
*2026-09-19*

Brüstle, *Derived-tame tree algebras*, Compositio Math. 129 (2001), Theorem 1.2:
a connected derived-tame tree algebra is determined **up to derived equivalence**
by three numbers — the vertex count, the corank of its Euler form, and the Dynkin
type of that form. Every LNA is a tree algebra with no extra hypothesis, and
Theorem 1.1 says derived-tame is exactly "Euler form positive semidefinite". So
for every LNA whose Euler form is non-negative, the class is decided with **no
search and no look at the relations at all**: invert the Cartan matrix, symmetrise,
take the rank and the Smith normal form.

**Implemented twice, independently, and the two agree.** `G = C^{-1} + C^{-T}`;
derived tame iff `G` is positive semidefinite; corank `= n - rank G`; Dynkin type
read off `(rank, product of the nonzero elementary divisors)` against the Cartan
determinants `A_r → r+1`, `D_r → 4`, `E_6 → 3`, `E_7 → 2`, `E_8 → 1`.

| | tame / total | classes | Coxeter polynomials per class |
|---|---|---|---|
| `n = 9` | 697 / 1430 | 6 | 1 each |
| `n = 10` | 810 / 4862 | 5 | 1 each |
| `n = 11` | 1514 / 16796 | 4 | 1 each |

At `n = 9` the six are `A_9` (128 members), `D_9` (145), `D̃_8` (38), `Ẽ_8` (377),
the tubular `(2,4,4)` (8) and `3033030` alone (corank 2, type `D_7`). That is
**exactly F-011's partition of the tame part**, arrived at without any mutation.

**It merges, and it never contradicts.** Against the full move orbits:

| | our orbits it merges | our orbits straddling two of its classes |
|---|---|---|
| `n = 9` | 2 into 1 | **0** |
| `n = 10` | 2 into 1 | **0** |
| `n = 11` | 2 into 1 | **0** |

One new merge at each length, every time inside the corank-0 type `D_n` class, and
never a contradiction — which is the check that matters, since a single orbit
split across two of Brüstle's classes would falsify either the theorem as stated
here or our moves.

**Why this is the first of its kind here.** Everything else we have is one-sided.
The quipu theorem names a class but only under almost separate relations; the
Coxeter polynomial rules a merge out and never in; the moves prove equivalence and
never inequivalence; the piecewise-hereditary criteria certify exclusion only.
Brüstle's is **two-sided on its domain**: two derived-tame LNAs agreeing on all
three numbers *are* derived equivalent, and disagreeing on any one they are *not*.

**Where it stops, and it stops hard.** The derived-**wild** LNAs, which are the
majority and are growing: 733 of 1430 at `n = 9`, 4052 of 4862 at `n = 10`, 15282
of 16796 at `n = 11`. Theorem 1.2 says nothing about any of them, and they are
where the open problem lives — F-046's three merges are all in there, and so is
the cospectral pair of F-010. What Brüstle does do at the boundary is split those
wild algebras off from the tame class that shares their Coxeter polynomial: at
`n = 10`, `T^10 + T^9 + T + 1` is carried by the 321-member tame `D_10` class and
by `34504030` and `50505000`, and the Euler form separates the first from the
other two. It says nothing about how those two stand to each other — F-037 has
that, and they are one class.

Note "non-negative" means **positive semidefinite**, not weakly non-negative; a
search over the positive cone does not see these.

Reproduce: `bruestle.py` and `bruestle_vs.py` in the sweep's scratchpad. E-040,
`research/literature/bruestle-derived-tame-tree-algebras.md`.

---

## F-044 — Happel–Seidel's table and our own naming agree, on both halves
*2026-09-19*

`research/literature/happel-seidel-piecewise-hereditary-nakayama.md` had to be
written from **secondary sources** — the paper is journal-only and was not
reachable — so the table in it is at second hand, and a table at second hand that
nothing checks is a liability. It checks out, on both halves, by two independent
routes the repo already had.

The setting is `N(n, r) = kA_n / rad^r`, our `[r] * (n - r) + [0] * (r - 2)`.

**Module type — the rows the table says are derived equivalent to a hereditary
star `T(a,b,c)`.** Compare `invariants.coxeterKey` of the LNA against
`treeSearch.treeCoxeterKey` of the star:

| | claimed | |
|---|---|---|
| `N_5(3)`, `N_6(4)`, `N_7(5)`, `N_10(8)` | `[2,3,r-1]` | match |
| `N_6(3)`, `N_9(5)` | `[2,3,r]` | match |
| `N_7(3)`, `N_8(3)`, `N_8(4)`, `N_9(3)`, `N_10(5)` | as tabulated | match |

**11 of 11 match.** A twelfth row, `N_11(6)` against `T(2,3,8)`, does **not**
match — and it is not one the table claims; it was invented here as a control, by
extending the `N_{r+2}` and `N_{r+3}` families to an `N_{r+5}` that the table does
not have. The table is not over-claiming, and the check has teeth.

**Sheaf type — the rows the table says are piecewise hereditary but derived
equivalent to a *canonical* algebra, so in no quipu class.** Compare the LNA's
Coxeter polynomial against `piecewiseHereditary.canonicalWeightType`:

| | Happel–Seidel | ours |
|---|---|---|
| `N_9(3)`, `N_9(5)`, `N_9(6)`, `N_9(7)` | `C(2,3,5)` | `(2,3,5)` |
| `N_9(4)` | `C(2,4,4)` | `(2,4,4)` |
| `N_10(3)`, `N_10(6)` | `C(2,3,6)` | `(2,3,6)` |
| `N_10(4)` | `C(2,4,5)` | `(2,4,5)` |
| `N_11(3)`, `N_11(6)`, `N_11(7)` | `C(2,3,7)` | `(2,3,7)` |
| `N_8(4)` | `C(2,3,4)` | `(2,3,4)` |

**12 of 12 agree.** And the split between the two halves falls exactly where the
theory says it must: enumerating every tree of the order against each polynomial,
a tree exists precisely for the **domestic** weight types (`1/p + 1/q + 1/s > 1`:
`(2,3,4)`, `(2,3,5)`) and for none of the tubular or wild ones (`(2,4,4)`,
`(2,3,6)` tubular; `(2,4,5)`, `(2,3,7)` wild). That is the trichotomy
`piecewiseHereditary` is built on, arrived at from the other end.

**What this is worth.** Two things, and neither is a new merge.

*The summary can be trusted.* Twelve published weight types and eleven published
tree types, none of which came from this codebase, reproduced by it.

*`canonicalWeightType` can be trusted.* It is the step of the pipeline with the
least independent support — a positive identification made from a Coxeter
polynomial, which is not a complete invariant — and this is the first time its
output has been checked against values published by someone else. Twelve for
twelve.

**What it is not.** No row here is a class the pipeline could not already name;
`N(n,r)` is a radical power, and the classification names those already. The
value is the audit, not the coverage. And the agreement is on Coxeter
polynomials, so it inherits their weakness: it confirms that we and Happel–Seidel
compute the same invariant and read it the same way, not that either reading is a
proof.

Reproduce: the two comparisons are in E-040.
`research/literature/happel-seidel-piecewise-hereditary-nakayama.md`.

---

## F-043 — The literature's radical-power merge is already in the move table
*2026-09-19*

arXiv:2302.02880 (Ueda) proves a triangle equivalence `per N(n, l+1) -> per N(n, l)`
for the radical-power Nakayama algebras `N(n, l) = kA_n / rad^l`, whenever
`n = p(p+1)q + p(p-1)r` and `l = (p+1)q + pr` for integers `p >= 2`, `q >= 1` and
`r >= 0` (or `p = 2` and `r` a half-integer). In our notation `N(n, l)` is the LNA
`[l] * (n - l) + [0] * (l - 2)`, so this is a statement about a thin but infinite
family of the rows we classify, and it is the only merge *between different
relation lengths* that any paper hands us outright.

**Every instance of it at `n <= 16` is joined by the moves we already have**, and
without needing the free move: 15 parameter triples fall in range, and for each,
`freeMoves.movesJoin(n, N(n,l), N(n,l+1), free = False)` finds a row both orbits
reach, checked from both ends with `orbitOf(..., target = meet)`.

**Which move does it matters, and the answer is not the rule table.** Asked with
each move set in turn, at every instance from `n = 10` up:

| | `n = 10` to `16` |
|---|---|
| the rule table alone | **no join, at any of them** |
| table + `edgeMoves` | no join |
| table + edges + `doubleMutation` | **join, at every one** |

So it is `proposition:doubleMutation` of arXiv:2310.08346 (F-032) that carries
these, and the table reaches them only below `n = 10`. Anyone re-measuring this
with `lnaMoves.closureUnderMoves` and the default rules will get "not joined" from
`n = 10` on and be right; the claim above is about the whole move set.

| | `N(n, l)` | `N(n, l+1)` | they meet at |
|---|---|---|---|
| `n = 6` | `3330` | `4400` | `4030` |
| `n = 8` | `555000` | `660000` | `605000` |
| `n = 12` | `4444444400` | `5555555000` | `4444555000` |
| `n = 16` | `10,10,10,10,10,10,0,…` | `11,11,11,11,11,0,…` | `80700558000000` |

(the full run is E-039; `10,10,…` is a comma list because a relation of ten arrows
has no single digit.)

**Two things this is worth, and one it is not.**

*It is an independent check on the move table.* The table was discovered by
search and verified against the Coxeter polynomial (R-005); here a published
theorem, proved by tilting objects and exceptional sequences rather than by
mutation, predicts 15 specific merges and the table produces all 15. Nothing in
the derivation of the table knew about Ueda's paper.

*It says where the literature currently sits relative to us.* The most on-target
merge result found in the sweep is **subsumed** at every length we can compute
at. Radical-power LNAs are a one-parameter family; the rows that are actually
open are the ones with relations of several different lengths overlapping, and
no paper found so far speaks about those.

*It is not evidence that the table is complete.* `N(n, l)` has all its relations
the same length and packed against the source, which is the easiest shape for the
moves; F-042 shows that where a cluster sits is what decides reachability, and
these sit where reachability is easiest.

Reproduce: `ueda2.py` in the merge session's scratchpad, or `movesJoin` on any row
of the table. E-039, `research/literature/2302.02880-ueda-derived-equivalences-nakayama.md`.

---

## F-042 — Where a cluster sits decides it, and no bound on the overlap does
*2026-09-18*

*Renumbered at merge from `F-039`, which was taken on `main` first by an unrelated entry while this branch was open. Session logs and commit messages from the branch use the old identifier.*

The quipu theorem names the LNAs whose relations are *almost separate* -- every
consecutive pair sharing at most one arrow -- so the obvious generalisation to
look for is a weaker bound on the overlap: at most two shared arrows, or at most
so many pairs that share more than one. Crossing those two coordinates against
membership of a quipu class (the theorem's own rows, closed under the move table,
the edge moves and the double mutation) says the coordinates are wrong.

| (max overlap, how many overlaps `>= 2`) | `n = 9` | `n = 10` | `n = 11` |
|---|---|---|---|
| (0, 0) and (1, 0) -- almost separate | 0 / 610 | 0 / 1597 | 0 / 4181 |
| **(2, 1)** | **1** / 300 | **12** / 954 | **84** / 2939 |
| (2, 2) | 0 / 132 | 20 / 483 | 169 / 1671 |
| (3, 1) | 0 / 90 | 21 / 300 | 152 / 954 |
| (4, 1) | 0 / 25 | 4 / 90 | 51 / 300 |

(outside a quipu class / total). The cell `(2, 1)` is the smallest possible step
past the theorem -- one pair of relations sharing two arrows, everything else
separate -- and it already contains outsiders at every length. No bound on the
overlap, however generous, cuts them off, and no bound on how many overlaps
exceed one does either.

**What does decide it.** The outsiders with the most free arrows are the same two
tiny configurations at every length from 10 to 12: a four-arrow relation followed
by a five-arrow one sharing three arrows (`45`), and `504`. Placing the `45` core
at every offset:

| | gap to the source: 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| `n = 9` | inside | inside | inside | | | | |
| `n = 10` | inside | **outside** | inside | inside | | | |
| `n = 11` | inside | **outside** | **outside** | inside | inside | | |
| `n = 12` | inside | **outside** | **outside** | **outside** | inside | inside | |
| `n = 13` | inside | **outside** | **outside** | **outside** | **outside** | inside | inside |

Read along a row: the moves carry the core to an almost separate LNA exactly when
it sits **against the source, or within one arrow of the sink**, and nowhere else.
The outside band grows one place longer with every vertex added, and the offsets
1 to 5 at `n = 14` are outside too. The overlap is three at every one of those
placements, the relations number two, and nothing about the cluster changes --
only where it is.

**And it is this core, not the overlap.** The same slide for every other small
core, at `n = 12` and `n = 13`: `33`, `44`, `34`, `43`, `54`, `55` and `333` are
inside at **every** placement. `55` overlaps in four arrows and is always inside;
`45` overlaps in three and is not. Even the direction matters -- `54` is inside
everywhere, `45` is not -- and that is not an inconsistency but the opposite
algebra at work: reversing the arrows sends the `45` core at offset `o` to the
`504` core at offset `n - 7 - o`, and `504` is outside at exactly the offsets that
sends them to. The outsiders with room are one family and its opposite, and the
size of the overlap does not pick them out.

**The short-quiver artefact, caught in the act.** At `n = 9` the core fits at
three offsets and **all three are inside**: `0450000` reaches an almost separate
LNA in 21 rows. Add one vertex at the far end and the same core, `04500000`, never
reaches one. So a length-9 census of what the moves reach overstates their reach,
and it does so for a configuration with two relations and a single overlap -- not
for an exotic one. F-040's warning has a concrete instance.

**The family.** `0^a 4 5 0^b` with `a >= 1` and `b >= 2` is outside the move
closure at every length tested, with as much free space on either side as wanted.
Anything that generalises the quipu theorem has to either place these or exclude
them by where they sit; a condition on the relations alone cannot see the
difference between `0450000` and `04500000`.

**What is measured.** Membership here is a *certificate*: reaching an almost
separate LNA proves the class is a quipu class. Failing to reach one is not a
proof that no derived equivalence exists -- it is a statement about this move set,
the same caveat `orbitOf` carries.

`freeMoves.movesFrom`, `overlap.overlapProfile`. E-037.

---

## F-041 — Meeting in the middle settles the free move at n = 8, and nearly at 9 and 10
*2026-09-17* · *re-measured 2026-09-19 under the guarded search, unchanged at `n = 8` and `n = 9` — E-038*

*Renumbered at merge from `F-038`, which was taken on `main` first by an unrelated entry while this branch was open. Session logs and commit messages from the branch use the old identifier.*

H-012 asks whether deleting a relation of two arrows -- a derived equivalence by
`corollary:lengthtworelations` -- is also a **mutation** equivalence. The sharper
question is one relation at a time, since the whole strip is a composition of
single deletions, and the sharper tool is to stop insisting that one algebra
reach the other.

**`search.meetingPoints`.** Every quiver a search passes through is in the class,
not just the lines it lands on. Two searches that arrive at the same quiver have
joined their algebras, at **twice the depth for the same cost**. Vertex labels do
not move under mutation, so two algebras on the same vertices meet on the nose
and the test is equality of labelled quivers, not isomorphism. Nothing in the
pipeline did this: `resolveMergeCandidates` collects only the lines a search
reaches and throws the rest of the tree away.

**What it settles.** Every single two-arrow deletion of every LNA:

| | `n = 8` | `n = 9` | `n = 10` |
|---|---|---|---|
| single deletions | 572 | 2002 | 7072 |
| joined by the known moves | 562 | 1937 | 6768 |
| joined by meeting in the middle, depth 3 | **10** | **52** | **206** |
| left open | **0** | 13 | 98 |

So **H-012 holds outright at `n = 8`**: every LNA there is mutation equivalent to
its stripped form, by composing single deletions. At `n = 9` thirteen deletions
are left, and eleven of those have both sides almost separate *with the same
quipu*, so the quipu theorem already calls them one derived class -- whether that
makes them one mutation class is a question about how that theorem is proved, not
about these searches. The two that are outside the theorem, `2302330 / 2300330`
and `3302302 / 3300302`, are the ones actually open, and they do not meet at
depth 4 either. Nor do any of the 98 left at `n = 10`: a depth-4 pass over them
joins none, so the extra depth buys nothing on either length and the pairs that
remain are either far further apart than 4 + 4 mutations or not joined at all.

**Why it works where a one-sided search does not.** The pairs it joins are joined
at 3 + 3 mutations: a one-sided search would need depth 6, where the pipeline
runs at 3 to 6 and the cost is exponential in the depth. Meeting in the middle
buys the same reach for the square root of the work.

`search.meetingPoints`, `search.quiversReachedFrom`. E-037.

---

## F-040 — What lengths 9 to 11 can and cannot show
*2026-09-17*

*Renumbered at merge from `F-037`, which was taken on `main` first by an unrelated entry while this branch was open. Session logs and commit messages from the branch use the old identifier.*

Every finding about the LNAs the quipu theorem misses rests on lengths 9, 10 and
11, and it is worth writing down exactly how narrow that evidence is. Counting
the **heavy clusters** of an LNA -- maximal runs of relations linked by overlaps
of two arrows or more, which is what puts an LNA outside the theorem:

| | `n = 8` | `n = 9` | `n = 10` | `n = 11` |
|---|---|---|---|---|
| LNAs with 0 heavy clusters | 233 | 610 | 1597 | 4181 |
| with 1 | 195 | 806 | 3148 | 11853 |
| with 2 | 1 | 14 | 117 | 761 |
| with 3 | 0 | 0 | 0 | 1 |
| **outside a quipu class** | 0 | 9 | 262 | 2647 |
| of those, with 2 clusters | -- | 0 | 2 | 42 |
| of those, with two clusters and a **free arrow between them** | -- | **0** | **0** | **0** |

So at every length the project has worked at, **an LNA outside a quipu class is a
single overlapping cluster**, and the few with two have them touching. Not one
has two clusters with a relation-free stretch between them. And the cluster is
usually against an end: at `n = 11`, 2119 of the 2647 have it touching the source
or the sink, 380 one arrow away, 123 two, 25 three.

**Where the missing configurations start.** Two heavy clusters with a free arrow
between them first fit at `n = 10` (5 LNAs, all of them in quipu classes). A
*barricade* -- two heavy clusters with a two-arrow relation walled in between
them, which is the shape H-012's doubts are about -- needs 4 + 1 + 2 + 1 + 4
arrows and so first fits at **`n = 13`**, two lengths beyond anything classified
here.

**What this does and does not undermine.** It does not touch F-033 or F-035, which
are statements about what was enumerated. It bears directly on anything phrased
as "every class", F-034 and H-014 above all: those say that every LNA outside a
quipu class reaches a quipu with relations, and every one of them is a single
cluster near an end. Whether that survives several clusters far from both ends is
untested and untestable at these lengths -- and the interaction between separated
clusters is exactly what a classification would have to handle in general.

`overlap.overlapRuns`. E-037.

---

## F-039 — Most of F-038's bad steps were bad *measurements*, and naming the arrows fixes them
*2026-09-18*

F-038 walked every LNA's search tree comparing `coxeterKey` at each node against
the start's, found nodes where it had moved, and split them into "parallel
arrows, a limitation of the model and harmless" and "clean, and the real fault".
**Both halves were misread, and in the same way: the key was being computed
wrong, not moved by the mutation.** Two independent mis-counts, one fixed by
naming arrows and one by closing the commutativity relations properly.

### 1. A parallel pair is two paths, and the Cartan matrix counted one

The procedure of arXiv:2112.08129 **produces parallel arrows**. Step 1 adds a
composite `alpha beta: h -> j` for every `beta: h -> i` and `alpha: i -> j`, and
`h -> j` may be an arrow already; step 3 adds one arrow `i* -> k` per relation
`i ~~> k`, and two relations may share both ends. Neither is a degenerate case:
the first one shows up five mutations out of an LNA at `n = 6`.

A path was a sequence of vertices, so two parallel arrows gave one path and the
Cartan matrix entry read 1 where it is 2. The smallest instance is the Kronecker
quiver — two arrows `1 -> 2`, no relations — whose Coxeter polynomial is
`(lambda - 1)^2` and which the repo read as `A_2`.

`3030` at `n = 6`, mutated at `[1, 3, 4, 1, 4]`, is the case E-033 found first.
The fifth step produces two arrows `1 -> 6`, one the composite of `1 -> 4` and
`4 -> 6` and one that was there before:

    key (1, 1, -1, -2, -1, 1, 1)  ->  (1, 1, 0, -1, 0, 1, 1)   counting vertices
    key (1, 1, -1, -2, -1, 1, 1)  ->  (1, 1, -1, -2, -1, 1, 1) counting arrows

Every step is admissible, and **the mutation was a derived equivalence all
along**. What was also lost there is the relation: `5 -> 1 -> 4 -> 6 = 5 -> 1 -> 6`
runs through the mutated vertex on one side only, so afterwards its two paths use
the two *different* arrows `1 -> 6` — and as vertex sequences they read alike and
the relation collapsed to nothing.

### 2. The cheap path count was wrong in two further ways, both about relations

Independent of parallel arrows, and the reason the *other* column of F-038's
table is not what it says either. The Coxeter key is read off
`invariants.integerCartanMatrix`, which counted paths rather than taking the rank
of the ideal, because a search calls it at every node and the exact route costs
about three times as much.

**It did not close the commutativity relations.** `paths.numberOfPathsUpToRels`
identifies paths by applying each of a subset of the two-path relations once, in
every order, and comparing canonical forms. That is not the closure, and a path
that is zero only through a **chain** of identifications was counted as nonzero.
From `40030` at `n = 7` by `[1, 4, 2, 5, 2]` — acyclic, no parallel arrows, one
of the four nodes F-038 called "clean" at that depth:

    1,4,2,5,7  ~  1,6,2,5,7   by  1,4,2 = 1,6,2
               ~  1,6,2,7     by  6,2,5,7 = 6,2,7
               ~  1,4,2,7     by  1,4,2 = 1,6,2   = 0  by  4,2,7 = 0

so `dim e_7 A e_1 = 0` and the old count said 1.
`arrowPaths.homDimensionByClosure` takes the closure to a fixed point and agrees
with the exact answer here.

**And it has no reading of a relation with three or more paths at all.** A
one-path relation is "this is zero" and a two-path relation is "these two are the
same"; a *sum of three* is neither, and both the old count and the new closure
ignore it outright. **Step 4 of the procedure produces one at every vertex with
three arrows out**, so this is not exotic. From `34400` at `n = 7` by
`[1, 3, 4, 2, 2, 1]`:

    arrows      (1,2) (1,4) (1,6) (2,7) (3,1) (4,7) (5,1) (6,7)
    relations   -(1,2,7) + (1,4,7) + (1,6,7) = 0,   (3,1,6) = 0,   (5,1,2) = 0

The three paths `1 ~~> 7` span two dimensions, not three. The cheap count reads
3, the key moves from `(1, 1, 0, -2, -2, 0, 1, 1)` to
`(1, -1, -4, -8, -8, -4, -1, 1)`, and the rank over the ideal gives the starting
key back. **No amount of closure fixes this**, and it is the bulk of what F-038
counted as "clean" at `n = 7` depth 6.

**So the key is exact now, except where the cheap route is provably right.**
`arrowPaths.isMonomial` is the condition: with every relation a single path, a
path is zero exactly when it contains a generator and counting is the dimension.
Every LNA, every tree and every quipu with zero relations is of that kind, so the
seeds and the recorded answers still take the cheap route; a mutated quiver
generally is not, and takes the rank.

### 3. What the two fixes do to the sweep

Every LNA of the length, searched with `coxeterGuard = False` so the corrupt
region is visible, `coxeterKey` compared at every node:

| | | nodes | wrong key | parallel | clean | lines |
|---|---|---|---|---|---|---|
| `n = 6`, depth 5 | before | 14,693 | 4 | 4 | 0 | 1,789 |
| | after | 14,701 | **0** | 0 | 0 | 1,789 |
| `n = 7`, depth 5 | before | 94,446 | 79 | 75 | 4 | 7,175 |
| | after | 94,498 | **0** | 0 | 0 | 7,175 |
| `n = 7`, depth 6 | before | 336,760 | 339 | 277 | 62 | 20,683 |
| | after | 337,360 | **10** | 0 | 10 | 20,683 |

The node count *rises*, by 8, 52 and 600, because a parallel-arrow node is no
longer terminal and the search descends from it. **Not one line is lost or gained**
at any of the three, which is the acceptance test F-038's own fix was judged by.

The ten that survive at depth 6 are one bad step and its descendants: every one
is reached from `30330`, the relation dual of `33030`, along the `[4, 1, 3, 1, …]`
family of paths, and they are section 4.

(F-038's table counts 25,398 nodes and 3,263 lines at `n = 6` depth 5 against
14,693 and 1,789 here; that sweep searched each LNA *and its relation dual*, this
one searches each LNA. The two columns to compare are the wrong-key counts and
the before/after within each row.)

### 4. What still moves the key, and why the guard stays

**A real failure remains.** From the relation dual of `33030` at `n = 7` by
`[4, 1, 3, 1, 3, 3]` — F-038's own smallest clean case — the key moves at step 6
under the arrow model too, and the cheap and exact Cartan matrices agree there,
so it is the algebra that changed and not the measurement:

    before   arrows (1,4) (1,6) (2,5) (3,7) (4,3) (5,1) (6,3)
             relations  [1,4,3,7] = [1,6,3,7],  [2,5,1],  [5,1,4,3] = [5,1,6,3]
    after    arrows (1,4) (1,6) (2,5) (4,7) (5,1) (6,7) (7,3)
             relations  [1,4,7] = [1,6,7],  [2,5,1],  [4,7,3],  [6,7,3]

    key  (1, 1, -1, -1, -1, -1, 1, 1)  ->  (1, 1, 0, 0, 0, 0, 1, 1)

So R-012 stands and `coxeterGuard` stays on: the criterion still rules mutation
out rather than in, and a step it admits can still fail to be a derived
equivalence. What changes is **how much** of the corrupt region is real: none of
it at `n = 6` and `n = 7` to depth 5, and at `n = 7` depth 6 ten nodes out of 339,
all of them this one step and what follows it.

### 5. What was fixed, in the code

* `arrowPaths` — the model. An arrow is `(tail, head, key)`, a path is a tuple of
  arrows, the empty tuple is a trivial path, composition is concatenation. The
  ideal arithmetic is the same linear algebra as `relationAlgebra` over a
  different basis, and reuses its row reduction.
* `procedure` — steps 1 to 7 per *arrow*. Step 5 divides by the arrow `alpha`,
  not by its target. Step 7 reads each candidate's first arrow back as the
  relation it came from rather than skipping the target when two relations share
  it, and reads its tail back into the old quiver (a composite is the two old
  arrows) before testing membership in `I`. Step 6 and the carried-past relations
  name the composite arrow they use.
* `procedure.isMutable` — no longer refuses a quiver with a parallel pair. The
  refusal said in as many words that it was a restriction of the model and not of
  the procedure; there is nothing left to restrict.
  `allowParallelArrows = False` keeps the old gate for measurement.
* `invariants` — the Cartan matrix, exact and cheap, counts arrow paths.
* `search` — the illegal-relation check is over arrow relations.
  `paths.isIllegalRelation` read a commutativity relation between two parallel
  paths as a repeated path and discarded the mutation, so the branch died even
  before the gate refused it.
* `pathAlgebra` — `arrowRels` carries the arrow relations, and the dual carries
  them across `(tail, head, key) -> (head, tail, key)`, without which a
  parallel-arrow algebra could not be left-mutated at all.

`rels` stays the table's key and is a **lossy projection**: a relation between
two parallel paths projects to the same vertex sequence twice, and reading that
back as this repo reads a two-path relation gives `p - p = 0` — no relation at
all. `relationsFrom` checks the cache describes `rels` and names arrows the
quiver has before trusting it, and raises rather than guess when asked to lift a
vertex sequence along a parallel pair.

**What this does not settle.** Whether the region beyond a parallel-arrow node
*reaches* anything is open and is H-016: at `n = 6` and `n = 7` to depth 6 the
search reaches exactly the same LNAs with the gate allowing parallel arrows and
with it refusing them, which is a weak negative -- those lengths need no search
at all (F-021) and a return from the region costs more depth than was searched.
There is also no canonical form for a quiver with parallel arrows, so two such
algebras reached by different routes cannot be compared; nothing in the pipeline needs to today, because every answer is
recorded at a line and a line has no arrow to spare for a parallel pair. And
step 3's *cyclic* case is still not implemented (F-002) — the model can now state
its answer, which is a precondition, and the gate still refuses a loop.

Reproduce: `tests/test_parallel_arrows.py`, and the sweep is E-035.

---

## F-038 — The search performs mutations that are not derived equivalences, and walks on from them
*2026-09-18*

`mutationSearchDepthFirst` gated every step on `mutationIsPossibleAtVertex` and
nothing else. That gate **rules mutation out, not in** — the theorem's hypothesis
is `Hom(P_i*[1], Λ) = 0`, on the algebra, and arXiv:2112.08129 says plainly that
this is in general not equivalent to a condition on the quiver. R-005 recorded
exactly this for rule discovery, where 38 "rules" passed on admissibility alone
and their orbits had the wrong Coxeter polynomial 6561 times out of 8388, and
made `lnaMoves.verifyMove` require three things. The search only ever asked for
one of them.

**Measured.** Walking every LNA's search tree and comparing `coxeterKey` at every
node against the start's:

| | nodes visited | wrong key | parallel arrows | oriented cycle | **clean** | lines reported | lines wrong |
|---|---|---|---|---|---|---|---|
| `n = 6`, depth 5 | 25,398 | 4 | 4 | 0 | **0** | 3,263 | 0 |
| `n = 7`, depth 6 | 609,474 | 604 | 507 | 0 | **97** | 37,911 | 0 |
| `n = 8`, depth 5 | 1,093,976 | 1,204 | 1,030 | 0 | **174** | 55,175 | 0 |

**Two mechanisms, and only one of them is dangerous.**

*Parallel arrows* are a limitation of the model, not of the procedure: a path is
a sequence of vertices and cannot say which of two arrows it uses, so the Cartan
matrix is misread. These are harmless to answers. `procedure.isMutable` refuses
mutation at *any* vertex of a quiver that has parallel arrows anywhere, so such a
node is terminal; and a quiver on `n` vertices with `n - 1` arrows two of which
are parallel cannot have a path of length `n - 1`, so it can never be mistaken
for a line either.

*The clean ones are the real fault.* Acyclic, no parallel arrows, Coxeter key
moved — and since nothing about them refuses mutation, the search descends
straight through and everything below is outside the class. The smallest is at
`n = 7`, from the relation dual of `33030` by `[4, 1, 3, 1, 3, 3]`. Every step is
admissible and no relation is illegal. At step 6, mutating at vertex 3:

    before   arrows (1,4) (1,6) (2,5) (3,7) (4,3) (5,1) (6,3)
             relations  [1,4,3,7] = [1,6,3,7],  [2,5,1],  [5,1,4,3] = [5,1,6,3]
    after    arrows (1,4) (1,6) (2,5) (4,7) (5,1) (6,7) (7,3)
             relations  [1,4,7] = [1,6,7],  [2,5,1],  [4,7,3],  [6,7,3]

    key  (1, 1, -1, -1, -1, -1, 1, 1)  ->  (1, 1, 0, 0, 0, 0, 1, 1)

The commutativity relation `[5,1,4,3] = [5,1,6,3]` is simply gone.

**No reported answer was wrong at `n ≤ 8`** — across 96,349 lines collected in
the three sweeps above, every one carried the starting Coxeter key. **At
`n = 10` and depth 8 it does reach an answer**, and that is what the ALARM of
E-032 was. Searching from the relation dual of `03033030`, four paths of the form
`[4, 6, 4, 6, 9, 4, 4, 6]` report the line `30233330`, which is in a different
class. Steps 1 to 6 hold the key, **step 7 — a second consecutive mutation at
vertex 4 — moves it**, and step 8 lands back on a line. Every step is admissible
and every quiver acyclic with no parallel arrows.

That is why the small sweeps came up clean: the corruption has to be entered and
then *returned from* to a line, which takes more depth than `n ≤ 8` was searched
to. The failure is not rare at depth, and it is invisible without the invariant.

*Amended 2026-09-18 → F-039, R-013.* **The table above does not count what this
finding says it counts.** A node is in it because the Coxeter key *as computed*
differed from the start's, and the key was being computed wrong in two ways: a
parallel pair of arrows counted as one path, and the cheap path count did not
close the commutativity relations to a fixed point. Correcting both takes the
wrong-key counts at `n = 6` and `n = 7` to depth 5 to **zero** — every one of
them, the "clean" ones at that depth included, was a mis-measurement of a
mutation that was a derived equivalence. The claim in the title still stands:
the smallest clean case here, from the relation dual of `33030` by
`[4, 1, 3, 1, 3, 3]`, still moves the key under the corrected count, so the guard
is still needed. What is retracted is the count and the reading of the parallel
half as "a limitation of the model and harmless to answers" — it was a
limitation of the model that made the search *refuse correct mutations*.

**The fix, and what it costs.** `mutationSearchDepthFirst` now takes
`coxeterGuard`, on by default: a step whose `coxeterKey` differs from the start's
is not taken. It is R-005's third requirement applied per step.

* It removes the corruption entirely — the wrong-key counts above go to **zero**.
* It changes **no answer**: over every LNA at `n = 6` and `n = 7` to depth 5,
  searched from the member and from the dual, not one line is lost and not one is
  gained.
* It costs **1.85×** at both lengths.

Reproduce: the guard is `search.mutationSearchDepthFirst(..., coxeterGuard =
False)` for the old behaviour, and the sweeps are in E-033.

**A second bug, real but dormant.** `search.py` discarded a mutation yielding an
illegal relation with `break` where it meant `continue`, abandoning every
remaining vertex at that node — and since the loop runs over `reversed(vertices)`,
that is every lower-numbered one. It is now `continue`. It never fired in the
overnight run: `isIllegalRelation` prints when it triggers and all three logs
contain zero such lines over 121 core-hours, so no result of E-032 is weakened by
it.

E-033, E-034, R-012.

---

## F-037 — The pair the Coxeter polynomial cannot separate is one class, and there is a path to prove it
*2026-09-18*

At `n = 10`, `34504030` and `50505000` are the two LNAs whose Coxeter polynomial
is `T¹⁰ + T⁹ + T + 1`. That polynomial is also carried by a quipu class, which is
the point `remark:Coxeter` of arXiv:2310.08346 makes: a polynomial match is
necessary and not sufficient, and for these two the invariant can say nothing at
all. H-013 recorded them as singletons under every move known — neither the rule
table, the free move, the edge moves nor the double mutation joins them, and
neither is certified non-piecewise-hereditary by our criteria — and predicted no
link to depth 8.

**They are derived equivalent, by seven mutations.** A depth-7 search from
`34504030` collects 20 lines, **19 of which are `50505000`**, and the reverse
search finds it too (E-032). Five of the shortest paths, replayed one mutation at
a time:

    [4, 2, 1, 1, 2, 2, 4]
    [4, 1, 2, 2, 2, 4, 1]
    [4, 1, 2, 2, 2, 1, 4]
    [4, 1, 2, 2, 1, 2, 4]
    [4, 1, 1, 2, 2, 2, 4]

Each was checked at every intermediate quiver for four things: the mutation was
admissible by `mutationIsPossibleAtVertex`, no relation was illegal, **no
parallel arrow or oriented cycle appeared**, and the Coxeter key did not move.
All five pass on all four counts. The last two conditions are not ceremony —
F-038 is a failure of exactly that kind, and this verification is what
distinguishes a real link from one the search fabricated.

**What it settles.** The 12 leftover orbits at `n = 10` fall to at most 10
classes, so the derived classes at `n = 10` number between 43 and 46. More
importantly it answers a question H-013 framed as needing a *new invariant*: the
pair the polynomial cannot separate did not need separating. Before reaching for
τ-periodicity or Hochschild cohomology on a pair like this, search it.

**Where it does not reach.** Seven mutations is beyond every depth the pipeline
runs at by default (`classify.py` starts at 6 and decays), which is why an orbit
this small sat unresolved. And the other leftover group at `n = 10` — the four
orbits on `(λ-1)²(λ+1)²(λ²+λ+1)(λ⁴-λ²+1)`, every member certified not piecewise
hereditary — stayed apart through depth 8, so this is not a general collapse.

Reproduce: `python merges.py 10 --depths 5 6 7 8`, or the path check directly
from `34504030` at depth 7.

E-032, E-033.

---

## F-036 — Reorienting a relation-free tree is a sequence of mutations, and the hereditary form is a mutation invariant
*2026-09-17*

The classification names a class by the tree of any relation-free quiver its
search reaches, and merges two classes that reach the same tree. The
justification written down for that has always been the *derived* statement --
two tree algebras are derived equivalent exactly when the trees are isomorphic,
since the orientations are related by BGP reflections. But a class in this repo
is a **mutation** class, and two searches can reach the same tree in different
orientations. The step in between was never checked: **are those reflections
mutations?**

They are.

**Right mutation at a source is the reflection, exactly.** On a relation-free
quiver whose underlying graph is a tree, mutating at a source reverses precisely
the arrows at that vertex, creates no relations, and does not renumber. Left
mutation at a sink is its inverse. Every tree of orders 3 to 7 in every
orientation, at every source and every sink:

| | checked | not the reflection | refused by the procedure |
|---|---|---|---|
| right mutation at a source | 2339 | **0** | 0 |
| left mutation at a sink | 2339 | **0** | -- |

**Any two orientations are joined, and the sequence can be written down.** Take
an edge where they differ and a side of it holding no other differing edge; flip
every vertex of that side exactly once. Each edge inside the side is reversed
twice and comes back; the differing edge has one endpoint inside and is reversed
once. Ordering the flips along a topological order of the side makes every one of
them a source (or, in the mirror case, a sink), so every step is legal.
`reflections.reflectionSequence`.

Verified against the engine -- each step put to the procedure's own admissibility
test, right mutations directly and left ones on the opposite algebra, then the
result compared with the orientation asked for:

| order | reorientations verified | failures | sequence length, mean / max |
|---|---|---|---|
| 4-7 | 3840 (all pairs, sampled) | **0** | -- |
| 8 | 432 | **0** | 5.5 / 12 |
| 9 | 432 | **0** | 7.1 / 16 |
| 10 | 432 | **0** | 8.3 / 20 |

**Right mutations alone already suffice**, which matters because the search walks
only those: flipping sources, never sinks, still reaches every orientation of
every tree of orders 4 to 8. Flipping every vertex once in a topological order
reverses the whole quiver, and twice brings it back, so a source flip is
undone by source flips.

**But the distance is the point.** The worst case over all trees of an order:

| order | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|
| right mutations only | 4 | 6 | 9 | 12 | **16** |
| both directions | 3 | 3 | 6 | 6 | 10 |

A classification search runs at depth 6. At order 8 two orientations can be 16
right mutations apart, so a search cannot cross between them at any depth it can
afford -- and the construction supplies the path without searching at all.

**The merges really do rely on it.** Searching every LNA of a length to depth 4
and pairing up the ones that reach the same tree -- which is exactly what
`mergeReport` merges on:

| | `n = 7` | `n = 8` |
|---|---|---|
| LNAs reaching a relation-free quiver | 23 of 132 | 26 of 429 |
| trees reached | 4 | 5 |
| pairs reaching the same tree | 79 | 80 |
| pairs reaching **isomorphic quivers** | 61 | 65 |
| pairs whose orientations must be joined | **18** | **15** |

So in about a fifth of them the merge is a derived equivalence and nothing more
until the orientations are joined.

**And the merge is now constructive.** `reflections.mutationBridge` joins two
LNAs that reach a common tree by running one's path to the tree, reflecting onto
the other's orientation, and running the other's path backwards; the engine
confirms it lands on the second algebra up to relabelling.

| bridge, all at `n = 7` | to the tree | reflecting | back | total |
|---|---|---|---|---|
| `00030` to `33000` | 2 | **11** | 4 | 17 |
| `00300` to `03000` | 3 | **9** | 3 | 15 |
| `00400` to `40000` | 3 | **7** | 3 | 13 |
| `05000` to `50000` | 4 | 4 | 4 | 12 |

All four are pairs from the 18. The reflection legs are 4 to 11 mutations on
their own and the whole sequences 12 to 17, against a classification search that
runs at depth 6: these are merges that were being made already, correctly, and
whose path no search in the pipeline could have found.

`reflections`, `tests/test_reflections.py`. E-031.

---

## F-035 — A relation of two arrows is free on a line and on nothing else
*2026-09-17*

`corollary:lengthtworelations` of arXiv:2310.08346 says a relation of two arrows
does not change the derived equivalence type of a Nakayama algebra, and F-028
measured it: 4861 cases at `n = 3..10`, the Coxeter polynomial kept every time.
The natural guess is that the same holds wherever a two-arrow relation sits. It
does not, and the counterexample is tiny.

**`D_4` with one two-arrow relation.** Take `3 -> 2 -> 1` with `4 -> 2` and the
relation on `3 -> 2 -> 1`. Its Coxeter polynomial is `x^4 + x^3 + x^2 + x + 1`,
the polynomial of `kA_4`; deleting the relation gives `kD_4`, whose polynomial is
`x^4 + x^3 + x + 1`. So the deletion leaves the derived equivalence class, since
the polynomial is an invariant of it.

**How often, over the quipus.** Deleting every two-arrow relation from every
admissible monomial ideal on every orientation of every quipu:

| order | ideals with a two-arrow relation | polynomial kept | **changed** |
|---|---|---|---|
| 4 | 11 | 7 | **4** |
| 5 | 72 | 48 | **24** |
| 6 | 543 | 300 | **243** |
| 7 | 4160 | 2138 | **2022** |
| 8 | 34938 | 15337 | **19601** |

About half, and the half that keeps it is not an accident either: **on the
linearly oriented line the polynomial is kept in every case**, which is the
corollary itself and is the control this measurement needed.

**What it costs.** `quipuRelations`' `minArrows = 3`, which is what makes order
11 affordable, is a real restriction of the family and not a normalisation: the
ideals with a two-arrow relation are genuinely different algebras and a run that
leaves them out has not covered them. At orders 9 and 10 the runs use
`minArrows = 2` and the question does not arise.

`quipuRelations.freeRelationCheck`, `python families.py free 4 5 6 7 8`. E-030.

---

## F-034 — Quipus with relations carry every class the quipu theorem misses
*2026-09-17*

The quipu theorem is about quipu quivers with **no** relations. Relations are
what the mutation procedure spends its time on, and a walk from an LNA to its
quipu passes through quivers that have them, so there is no reason the shape
should stop mattering once it does. Enumerating every quipu of an order, every
orientation of it up to the tree's automorphisms, and every admissible monomial
ideal on the result, and matching the Coxeter polynomial against the LNAs of the
length:

| order | shortest relation | ideals walked | matching | polynomials of unclassified LNAs covered |
|---|---|---|---|---|
| 9 | 2 arrows | 370 483 | 3677 (2820 up to isomorphism) | **2 of 2** |
| 10 | 2 arrows | 3 411 263 | 306 624 | **7 of 7** |
| 11 | **3 arrows** | 5 465 194 | 1 246 011 | **20 of 20** |

Every Coxeter polynomial carried by an LNA that lies in no quipu class is also
carried by quipu algebras with relations, and by thousands of them: 1746 up to
isomorphism on `3033030` at `n = 9`, spread over 16 of the 18 quipu shapes of
that order; at `n = 11`, between 10 220 and 60 228 per polynomial, over 61 to 64
of the order's 64 shapes. The linearly oriented line is excluded throughout, so
none of this is an LNA matching itself.

The `n = 11` run allows only relations of **three arrows or more**, which is a
real restriction and not a normalisation (F-035), so its 20 of 20 is a statement
about that sub-family — and a stronger one for it, since the sub-family is
smaller.

**A polynomial match is not an equivalence, so the classes were walked as well.**
Every quiver a mutation search out of an LNA reaches is in its class by
construction, so the quipus among them are members with a path to prove it.

| | `n = 9` | `n = 10` | `n = 11` (sample of 200) |
|---|---|---|---|
| LNAs outside a quipu class | 9 | 262 | 2647 |
| **reaching a quipu with relations within 3 mutations** | **9** | **262** | **200 of 200** |
| reaching none | 0 | 0 | 0 |
| quipu algebras confirmed, per LNA: min / median / max | 8 / 16 / 18 | 3 / 20 / 84 | 5 / 27 / 117 |
| distinct quipu algebras confirmed | 178 (at depth 4) | 3510 | -- |

So it is not that *some* class outside the theorem has a quipu member: **every
LNA outside one has, and within three mutations** -- exhaustively at `n = 9` and
`n = 10`, and in a random sample of 200 of the 2647 rows at `n = 11`.

Every algebra reached this way is in the enumeration; the only quivers reached
and not enumerated were the linearly oriented lines, which are the LNAs
themselves — so the enumeration is complete in the sense that matters, and the
matches are not all coincidence.

**The smallest case is at order 4 and is not exotic.** `D_4` with one two-arrow
relation is `kA_4` after a single mutation (F-035). What is new is that the same
thing happens for the lines the theorem *cannot* name.

**What this does not say.** 178 of 2820 at `n = 9` are confirmed; the rest are
leads, and a lead is only a necessary condition met. Nor is there a *family* yet:
a theorem in the shape of `thm:QuipuToAn` would say which quipu-with-relations
goes with which LNA, and this says only that plenty of them exist. H-014.

`quipuRelations`, `python families.py quipus 9 --min-arrows 2`. E-030.

---

## F-033 — No tree outside the quipu shape carries an unclassified LNA's polynomial
*2026-09-17*

F-031 enumerated the trees of maximum degree three that are not quipus, up to
order 12, and found none sharing a Coxeter polynomial with any LNA. This is the
same question of **every** tree, degree four and above included.

| order | trees | not quipus | sharing a polynomial with a quipu | **leads** |
|---|---|---|---|---|
| 9 | 47 | 29 | 3 | **0** |
| 10 | 106 | 70 | 0 | **0** |
| 11 | 235 | 171 | 7 | **0** |
| 12 | 551 | 424 | 15 | **0** |

694 non-quipu trees over the four orders, and not one of them shares its Coxeter
polynomial with an LNA that lies outside a quipu class.

**The 25 that do share a polynomial are refuted without a search.** In every one
of those cases the LNAs under that polynomial are all ones the moves place in a
quipu class, and the tree is *cospectral* with that quipu — the same polynomial,
a different tree. Two path algebras of trees are derived equivalent exactly when
the trees are isomorphic, so the tree is not derived equivalent to the quipu, and
so not to the LNA either. The Coxeter polynomial matching is the whole of the
coincidence.

**What it is and is not.** Four orders, not a proof for all of them. Taken with
F-031 it says the hereditary side of the classification is closed as far as
anything has been able to look: if there is a second family behind the classes
the theorem misses, its quivers **have relations** — which is F-034.

`treeSearch`, `python families.py trees 9 10 11 12`. E-030.

---

## F-032 — The double mutation of arXiv:2310.08346 is the mechanism the rule table was approximating
*2026-09-17*

`proposition:doubleMutation`: for a relation `r: s → t` with a relation starting
at `s - 1` and none at `t - 1`, two left mutations at `t` add `s+1 → t+1` (unless
`t = n`), move the start of every relation starting inside `r` up one, and move
the end of every relation ending inside `r` up one. Stated in full, with the
proof's structure, in `literature/2310.08346-*.md`; implemented as
`quivermutation/doubleMutation.py`, the dual `R_s` taken through F-026.

**It holds against the engine everywhere it applies.** Every LNA of lengths 5 to
10, every relation satisfying the hypotheses, both directions, checked the three
ways `verifyMove` checks a rule — predicted LNA, every mutation admissible,
Coxeter polynomial fixed:

| n | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|
| confirmed | 18 | 68 | 250 | 922 | 3430 | 12868 |
| failures | 0 | 0 | 0 | 0 | 0 | 0 |

17556 in all. Included are the cases at `s = 1` with no companion, which the
paper does not state (its hypothesis is used only to keep the quiver a line, and
at the source it is one anyway): 1429 per direction at `n = 10`, no failures.
`tests/test_double_mutation.py` pins the counts.

**It is why no rule could state it.** Every relation it changes crosses `r`, so
every such relation straddles any window, and `matchesAt` refuses those. The
table holds only its instances with nothing crossing: `pairSlideRules` (588 of
588 applications at `n = 5..10` are one double mutation), `endPairCollapseRules`
(600 of 600), and `edgeMoves.sourceDoubling` (every one, `n = 7..9`). Of all
1844 table rules' applications at `n = 9`, 1602 of 3754 are a single double
mutation.

**What it reaches, with no search** (seeding from the quipu theorem, then orbits;
`python overlaps.py 9 10 --free --doubles --no-rules`):

| n | rule table + edges (mutation) | double mutation alone (mutation) | + free move (derived) | left | orbits left |
|---|---|---|---|---|---|
| 8 | 428 / 429 | **429 / 429** | 429 | 0 | 0 |
| 9 | 1292 / 1430 | 1397 | **1421** | 9 | 2 |
| 10 | 63 % | | **4600 / 4862 (94.6 %)** | 262 | 16 |
| 11 | 47 % | | **14149 / 16796 (84.2 %)** | 2647 | 86 |
| 12 | | | 42836 / 58786 (72.9 %) | 15950 | |

The rule table on top of the double mutation and the free move adds **nothing**
at `n = 10`: the same 262 rows in the same 16 orbits. The whole run takes seconds
where the table took minutes.

**At n = 9, 10 and 11 what is left is provably not a quipu class**, so seeding,
the free move and the double mutation together place every LNA that *is* in one:

- `n = 9`: the 9 rows are two orbits, exactly `3345000` (8 members, `C(2,4,4)`)
  and `3033030` (1, not piecewise hereditary) — F-011's two non-quipu classes,
  with nothing searched.
- `n = 10`: 16 orbits under 7 Coxeter polynomials. Six polynomials belong to no
  quipu of order 10. The seventh, `T^10 + T^9 + T + 1`, is carried by two
  singleton orbits `34504030` and `50505000`, which are precisely `Λ` and `Λ'` of
  the paper's `example:A10double`, proved not piecewise hereditary there by
  `τ^9(P_2) = P_2[1]`; neither of our criteria certifies them.
- `n = 11`: 86 orbits under 20 polynomials. One polynomial is shared with a
  quipu, and both orbits carrying it contain a certified non-piecewise-hereditary
  member (6 of 8 each), so neither is a quipu class.

**Caveats.** The free move is a derived equivalence only, so the partition with
it is of derived classes, not mutation classes (H-012). Orbits left over that
share a polynomial may still be one class. Closed under the relation dual as
well, which is free, `n = 10` has 12 orbits in 7 polynomial groups and `n = 11`
has 54 in 20, so 43 to 48 and 84 to 118 derived classes (H-013). The interior half of the move alone reaches much
less (273 / 429 at `n = 8`, 770 / 1430 at `n = 9`), so the ends are still where
the work is done (H-010, H-011). `E-029`.

---

## F-031 — The smallest tree that is not a quipu has ten vertices, and no LNA reaches it
*2026-09-16*

A quipu is a tree of maximum degree three whose degree-three vertices all lie on
one path. The natural question is whether the quipu theorem's shape is forced:
is some tree of maximum degree three that is *not* a quipu still derived
equivalent to a linear Nakayama algebra?

**Where the first one is.** Enumerating every tree of maximum degree three by
order, and testing each with `quipuForms.isQuipuByDegrees`:

| order | trees of max degree 3 | not quipus |
|---|---|---|
| 4–9 | 2, 2, 4, 6, 11, 18 | **0** |
| 10 | 37 | **1** |
| 11 | 66 | 2 |
| 12 | 135 | 8 |

So ten vertices is where they start, and at ten there is exactly one: the centre
`c` with three neighbours, each of which carries two leaves — `1 + 3 + 6 = 10`.
Its four degree-three vertices are `c` and the three neighbours, and a path
through the tree can hold `c` and at most two of the others, so no path holds all
four.

**No LNA of that length is derived equivalent to it.** The Coxeter polynomial is a
derived invariant, and derived equivalence forces the same number of simples, so
the candidates are exactly the LNAs of the same length. Comparing against all of
them:

| order | non-quipu trees | LNAs compared | sharing a Coxeter polynomial |
|---|---|---|---|
| 10 | 1 | 4862 | **0** |
| 11 | 2 | 16796 | **0** |
| 12 | 8 | 58786 | **0** |

Eleven trees, 80444 LNAs, not one match. The Coxeter polynomial is only a
necessary condition, so a match would have needed following up; its absence is
decisive on its own.

**What this is and is not.** It is not a proof that no non-quipu tree is ever
reached — that is a statement about all orders, and this is three of them. It is
evidence that the quipu shape in `thm:QuipuToAn` is not an artefact of how the
theorem was proved, and it closes the cheapest way the classification could have
been incomplete. The next orders are the obvious extension and cost only time.

E-027.

---

## F-030 — An unequally overlapping pair becomes a triple, for every pair of lengths
*2026-09-16*

Opening a relation into its square and walking the side it opens (E-027) turns
up one family in the **interior**, and it is parameterised by *both* relation
lengths — which is what H-008 asked for and what F-020's three families, all
about a two-arrow relation sliding, did not give.

    (0:shorter) (1:longer)  ->  (0:shorter + 1) (1:longer) (2:longer)

under two left mutations at the window's last vertex, for every `longer >
shorter >= 2`. The window is `longer + 2` arrows. Verified at lengths 8 to 13,
**no failures**: 1991 confirmations for `(2,3)`, 603 apiece for `(2,4)` and
`(3,4)`, fewer for the wider windows simply because a wide window has fewer
places to sit.

**It is floating**, so it holds in the interior as well as at an end — and it is
the first family found there that changes the *number* of relations. Read
backwards it takes a run of three relations to a run of two, which is the
mechanism behind F-022's observation that a run of three dissolves where an
isolated pair does not.

**It does not lower the overlap, and could not.** `(3,4)` has a maximum overlap
of 2 and `(4,4,4)` of 3, so the rewrite runs *up* the overlap scale. That is
consistent with H-010 rather than against it: the family says an interior pair
has somewhere to go, not that it can get anywhere useful.

**H-008's shape once more.** The table already contained `longer = 3` and
`longer = 4` and none of the rest, because `longer = 5` needs a window of seven
arrows and discovery never had one. Generating the family instead of listing it
takes the floating table from 364 rules to 414.

**What it does not buy.** Nothing, yet, for the coverage: at `n <= 10` the new
members are too wide to fit anywhere the narrow ones did not already reach, and
the orbit counts do not move. The value is that the family is now known in
closed form, and it will matter at the lengths where its members fit.

`pairToTripleRules`. E-027.

---

## F-029 — A relation at an end doubles, and the rule encoding cannot say so
*2026-09-16*

    A relation of `l` arrows at the source, with nothing starting at the
    vertices 2 .. l, gains a second relation of `l` arrows at vertex 2, under
    two left mutations at vertex `l + 1`.

Verified against the mutation engine over every LNA carrying the pattern,
checked all three ways `verifyMove` checks a rule — the predicted LNA, every
mutation admissible, the Coxeter polynomial held fixed. For the source doubling:
`l = 2..7` at lengths 8, 9, 10 and `l = 2..9` at lengths 11 and 12 (`l = 8` and
`l = 9` need a longer quiver than 10 to appear at all), **9855 confirmations and
no failures**, with no case where the predicted LNA failed to be admissible. The
dual at the sink, by F-026's relation dual, holds identically, and all four moves
— both doublings and both collapses — fire exactly **907 times apiece** over
lengths 7 to 10 with no failures, the equality between a move and its dual being
itself a check.

**Why no discovery run could have found it.** `describeLink` builds the window
from the relations that change and the vertices mutated, here the arrows
`1 .. l + 1`, and then refuses to describe the rewrite unless every relation
meeting that window lies inside it. A companion relation starting on arrow
`l + 1` and running out the far side breaks that condition, and the rewrite is
discarded. So this family is invisible to the encoding, not merely to the search
bound — which is why 1794 verified rules do not contain it, and why in `A_12`
the whole table offers *no* rewrite at all for `3003000000` while this gives one.

The strictness is not a mistake: it is what stops a rewrite depending on
something it does not describe, and R-009 is what happens without it. The right
fix is to let a description say that relations may cross the window's far edge
untouched. That is a change to `matchesAt`, `applyAt` and `describeLink`, and it
is NOTES backlog 27; until then the family lives in `edgeMoves` as explicit
rewrites rather than as table rows.

**What it is, in one line.** The doubling is the inverse of `endPairCollapse`,
released from needing the rest of the window to be empty. Its collapse direction
is *not* the naive inverse of its sequence: two left mutations at vertex `l + 1`
are undone by two **right** mutations at the **source**, because the procedure
relabels. Inverting the condition by hand instead of defining the collapse as
"the LNA whose doubling is this one" made a first version fire on a third of the
cases it should have, and that is how it was caught.

**What it buys.** With the edge moves alone, `A_9` goes from 380 move orbits to
237 and from 222 rows needing a search to 138; `A_10` from 2149 to 1607 and from
1807 to 1426. With F-028's free move as well, **`A_8` needs no search at all** —
21 orbits, nothing left over — and `A_9` falls to 77 orbits and 37 rows.

`edgeMoves`, `tests/test_edge_moves.py`. E-027.

---

## F-028 — A relation of two arrows is free, and that is worth more than the whole rule table
*2026-09-16*

`corollary:lengthtworelations` of arXiv:2310.08346 says a relation of **two
arrows** does not affect an algebra's derived equivalence type. It has been in
`research/literature/` since that paper was read and has never been used. It is
the single largest reduction available to this investigation, and using it costs
nothing.

**It is not a mutation.** Every rule in `lnaMoves` is a rewrite the engine has
been made to perform, with a sequence behind it. This has no sequence. It is a
statement about *derived* equivalence only, and may never be used to claim two
algebras are in the same mutation class. That is why it lives in `freeMoves`
and not in the move table.

### Two independent checks, since the corollary is quoted and not proved here

| check | cases | disagreements |
|---|---|---|
| Coxeter polynomial kept under stripping, `n = 3..10` | 4861 | **0** |
| quipu name kept under stripping, `n = 3..13` | 44320 | **0** |

The Coxeter polynomial is a derived invariant, so it had to be blind to these
relations; the quipu theorem of arXiv:2305.06642 is what the classification is
seeded from, and it is blind to them too. Neither is a proof — the first is a
necessary condition and the second is a theorem about a subclass — but between
them they would have caught a misreading of the corollary, and did not. The
paper's own classification table, transcribed in `tests/paper_classification.py`,
already writes only relations of length ≥ 3 and says why.

### What follows, and the third one is the useful one

**A relation of two arrows never overlaps a neighbour in more than one arrow.**
Admissibility makes starts and ends both strictly increase, so for `(s, 2)` and a
later `(t, l)` we have `s < t` and the overlap `s + 2 - t <= 1`; for an earlier
`(r, k)` we have `r + k < s + 2` and the overlap `r + k - s <= 1`. So a two-arrow
relation can never be part of what F-021 says a search still has to find.

**Therefore stripping never names anything by itself.** If an LNA has an overlap
of two or more it is between two relations of three or more arrows, both of which
survive stripping and stay consecutive. So `isAlmostSeparate` is exactly
preserved — measured, 0 exceptions at `n = 5..10` — and the count of LNAs that
stripping alone hands to the quipu theorem is **0 at every length**. Everything
the move buys is in *bridging*: an LNA no rule reaches strips to one the rules do.

**The reduced space is one vertex smaller.** Shortening every relation by one
arrow is a bijection from the stripped LNAs of length `n` onto *all* the LNAs of
length `n - 1` — a reduced LNA has every relation at three arrows or more, so one
comes off each and the starts and ends still increase, and its last entry is
always 0 because a relation starting at `n - 1` could only have two arrows.
Checked as a bijection, not just a count, for `n = 4..13`.

| n | LNAs | reduced | ratio |
|---|---|---|---|
| 8 | 429 | 132 | 3.3 |
| 10 | 4862 | 1430 | 3.4 |
| 12 | 58786 | 16796 | 3.5 |
| 14 | 742900 | 208012 | 3.6 |

So quotienting by the free move divides the whole space by the Catalan ratio,
tending to 4.

### What it does to the coverage, which is the point

Orbits under the 1794 rules, and under the rules plus the free move:

| n | LNAs | orbits, rules | orbits, + free | needing a search, rules | + free |
|---|---|---|---|---|---|
| 7 | 132 | 20 | **10** | 0 | 0 |
| 8 | 429 | 69 | **22** | 10 | **1** |
| 9 | 1430 | 380 | **91** | 222 | **53** |
| 10 | 4862 | 2149 | **746** | 1807 | **971** |
| 11 | 16796 | 9254 | **3734** | 8863 | **6262** |
| 12 | 58786 | 34392 | **14340** | 37366 | **29312** |

One line of code merges 20052 pairs at `n = 12` that 1794 verified rules — the
product of every discovery run this project has made, and the table as it stood
when this was measured, before F-030 added 50 more — do not. At `n = 9` it cuts
what a search must still place by 76%, at `n = 8` by 90%.

### Why the rules missed it, and how far they get on their own

The table can slide a two-arrow relation anywhere (F-020's lone slide) and can
delete one at an end **when it is alone there**: in `A_8`, `(0,0,0,0,0,2)` reaches
the empty relation set, `(0,0,0,2,0,0)` is only ever slid, and `(3,0,0,0,2,0)` is
not stripped even though its two-arrow relation can be slid to the sink. Against
an end the deletion is one mutation — E-027 finds `(L, …, 2) → (L)` by a single
left mutation at the sink, for every `L` it tried. So the free move is not exotic;
it is the end behaviour, released from needing the end.

**This does not make the free move redundant.** Whether every instance of it is
*also* a mutation equivalence is a separate and open question: sliding a two-arrow
relation to an end needs room and company it may not have. What is settled is that
the derived classes merge, and derived classes are what is being classified.

`freeMoves`, `overlaps.py --free`, `tests/test_free_moves.py`. E-027.

---

## F-027 — Leaving the line, a mutation of an LNA goes to a square with a side of two
*2026-09-16*

A rule of two or more mutations passes through quivers that are not lines. They
are not arbitrary. Walking every multi-mutation rule in the table one step at a
time and classifying each intermediate:

| intermediate | count |
|---|---|
| commutative square, sides 2 and 2 | 450 |
| sides 2 and 3 | 351 |
| sides 2 and 4 | 297 |
| sides 2 and 5 | 216 |
| sides 2 and 6 | 31 |
| sides 2 and 7 | 9 |
| sides 2 and 8 | 10 |
| a line again | 336 |
| anything else | **5** |

1296 rules, 1705 intermediates. **Every square has a short side of exactly two**
-- 1364 of them, and not one with a short side of three. The five exceptions
have a branching vertex and no matching join, so they are not squares at all.

**What the square is.** Mutating at the source vertex of a relation of k arrows
takes the relation `v -> v+1 -> ... -> v+k` and replaces it with a commutative
square: a new two-arrow path from a vertex to the relation's end, commuting with
the k-1 arrows still on the line. Right mutation at the start does it one way
round and left mutation at the end the other:

```
A_10, one relation on the arrows 3..7, right mutation at 3
    relations  (2,4,3)  and  (4,3,8) = (4,5,6,7,8)
                              ^^^^^ two arrows          ^^^^^^^ four arrows
```

so the zero relation has become a *commutativity* relation between a side of two
and a side of k - 1.

**Why it matters, and it is a construction rather than an observation.** The
rules were all found by search -- enumerate mutation sequences, describe what
recurs. The square says what the search is walking through: a rule is *open the
relation into a square, do something along its long side, close it back onto a
line*. That is the same shape as F-020's families, where the mutation count grows
one per arrow travelled, and it suggests building rules directly instead of
finding them: open at a chosen relation, walk the long side, and read off where
it closes.

**The caveat, from trying it.** A lone relation opened into a square closes only
by undoing itself -- searching four mutations from the 2-by-4 square of
`00500000` in A_10 finds nothing but `[-3]` back to where it started. The square
has to have something to interact with, which is the same lesson as F-023's
spectators: the interesting rules are the ones with a second relation in the
window. Constructing rules this way therefore means opening a square *and*
choosing the companion, which is a smaller search than the one being run now but
not a formula.

E-026.

---

## F-026 — A rule's dual is a rule, and the table was missing 410 of them
*2026-09-16*

The relation dual -- reverse every arrow of the line and renumber -- preserves
the derived equivalence class of any LNA, and left mutation at a vertex is right
mutation at that vertex of the dual. So a rule must carry over to the dual
picture, and the transform is mechanical.

**The transform.** For a rewrite on a window of `width` arrows:

* a relation covering the arrows `s .. s+a-1` covers `width-s-a .. width-s-1`;
* the vertex at offset `o` becomes the one at `width-o+2`, and a **right**
  mutation there becomes a **left** one, and the other way about;
* the sequence keeps its order, the dual being applied step by step;
* an anchor to one end becomes an anchor to the other.

`lnaMoves.dualRule`, and it is an involution.

**It holds.** Of the 1384 rules the table then held, **none** was its own dual and
**410** had a dual that was not in the table. Verified where each fits, at up to
four lengths apiece: **410 hold, 0 fail, 0 never apply.** The table is generated
closed under the dual now -- 364 floating rules and 1430 anchored, 1794 in all.

**It is not the transform E-019 refuted**, and the difference is worth keeping
straight. That one tried to read a rule's **inverse** off its window by reversing
the sequence and negating: it worked for 12 of 96 rules and was abandoned. This
is a symmetry of the problem rather than a shortcut, and it works for all of
them.

**What it is worth, honestly: very little coverage.** Closing the table adds 410
rules and moves the count by **7 rows at n = 10 and none at n = 9**. The orbits
those duals join were already joined another way. Its value is elsewhere:

* it is free, and a table that is not closed under a symmetry of the problem is
  wrong to leave that way;
* it halves what a search has to look for -- discovery could plant patterns at
  one end only and dual the results, at half the cost of E-025's two hours;
* and it is the correction to F-025, which claimed an asymmetry between the two
  ends of the quiver on the strength of comparing a pattern with itself rather
  than with its dual (R-011).

E-026. Tests: `test_the_dual_of_a_rule_reverses_the_window_and_turns_the_mutations_round`,
`test_the_table_is_closed_under_the_dual`.

---

## F-025 — The two ends of the quiver are not the same end
*2026-09-16* · **RETRACTED 2026-09-16 → R-011**

**The asymmetry is the pattern's, not the quiver's.** The mirror of a rule is its
relation dual -- reverse the arrows *and* exchange right mutation for left -- and
under that transform the sink rule below holds perfectly well at the source, on
the dual pattern. What was compared with it was the same pattern at the other
end, which is a different configuration. The rules and the family below stand and
are in the table; the conclusion drawn from them does not. R-011, F-026.


`endPairCollapseRules` (F-022) collapses a pair of relations of **equal** length
against either end, and the two directions are mirror images, as the relation
dual says they must be. For a pair of **unequal** lengths that symmetry breaks,
and only one end works.

**The probe.** Each of the four long unequal pairs the residue of E-024 is made
of, planted flush against each end of A_13 or A_14 and mutated within three
vertices, at three mutations:

| pattern | at the source | at the sink |
|---|---|---|
| `(1:3) (2:6)` | 2 LNAs, **nothing lower** | overlap 2 → **0**, via `[7, 6, 7]` |
| `(1:3) (2:7)` | 2 LNAs, **nothing lower** | overlap 2 → **0**, via `[7, 6, 7]` |
| `(1:5) (2:6)` | 2 LNAs, **nothing lower** | overlap 4 → 3, via `[7, 7]` |
| `(1:5) (2:7)` | 2 LNAs, **nothing lower** | overlap 4 → 3, via `[7, 7]` |

At the source nothing moves at all. At the sink the **shorter** relation, which
is the one that starts first, loses an arrow — and where it had three arrows to
begin with, losing one takes the overlap to zero and the LNA into the quipu
theorem's reach outright.

**Two of them as rules, and the anchor is not decoration.**

```
window 7 arrows at the right end:  (0:3) (1:6)  ->  (0:2) (1:6)     via [2, 2]
window 8 arrows at the right end:  (1:3) (2:6)  ->  (0:2) (2:6)     via [3, 2, 3]
```

4 confirmations each over the three lengths their windows fit in, no failures.
The identical rewrites stated as floating rules give **4 confirmations and 4
failures** apiece: they are true against the sink and false elsewhere, checked
rather than assumed.

**Why the asymmetry is not surprising once stated.** The pair is `(n, l)` and
`(n+1, m)` with `l < m`, so the two relations start one vertex apart but end
`m - l + 1` arrows apart. Flush against the sink it is the *ends* that are
pinned, and the two relations end at different places, so the configuration
there is genuinely different from the one at the source, where it is the
*starts* that are pinned and they are one apart either way. The equal-length
pair is exactly the case where the two descriptions coincide, which is why
F-022's family is symmetric and this is not.

**What it opens.** The rule is stated here for the two pairs probed, not for a
family in `(l, m)`. Sweeping `l` and `m` is what E-025's discovery run is for,
and the family is the thing to look for in its output. E-025.

---

## F-024 — The interior is emptier than it looked, and the boundary was doing the work
*2026-09-16*

E-021 measured what a heavily overlapping pair can be turned into and found it
frozen. The measurement was right and the label on it was wrong: those probes
allowed mutations at every vertex of A_13, ends included, because a margin of 6
around arrows 5 to 8 of a 13-vertex quiver reaches both of them and even a margin
of 3 reaches vertex 2, which rewrites arrow 1. Re-run with the quiver lengthened
so that the ends are genuinely out of reach, the picture changes in size but not
in conclusion — and the difference is the finding.

**The same pattern, the same margin, with and without an end in reach.**

`(1:3) (2:3)`, margin 3, in A_13 at offset 4 (arrows 1 to 12 rewritable, so both
ends in reach) against A_21 at offset 8 (arrows 5 to 14 only):

| mutations | LNAs reached, ends in reach | LNAs reached, genuine interior |
|---|---|---|
| 3 | 8 | **2** |
| 4 | 14 | **4** |
| 5 | 22 | **4** |
| 6 | 36 | **6** |

In a genuine interior the pair is not so much frozen as nearly immobile: six
mutations reach six LNAs, and the fifth mutation buys nothing at all. Five of
every six LNAs the earlier probes reported were reached with the help of an end.

**And the overlap never comes down in the interior, now checked to six.** Every
one of those 2, 4, 4 and 6 still has two relations sharing two arrows, the
depth-6 run taking 2583 seconds to say so. That is H-010's claim tested two mutations
deeper than before, and it survives.

**With an end in reach, six mutations do pull the pair apart** — the one thing
that has ever done so:

```
A_13:  00003300000  ->  30000020000   via [-8, 5, 4, 3, 2, -6]
```

overlap 2 down to 0. The middle of that sequence is `5, 4, 3, 2`: four mutations
walking down the quiver, one vertex at a time, until the relation is at arrow 1
and there is no further to go. It is **not translation invariant** — shifted by
1, 2, 3 or 4 vertices in a quiver lengthened to match, it does not even produce
an LNA, let alone the shifted answer. So it is not a rule that happens to need
six mutations; it is the boundary, reached the long way round.

**What this settles and what it costs.** It settles that the escape route is the
one H-011 names -- walk the run to an end -- and that it is a single mutation
sequence, not only a composition of table rules. It costs E-021 its framing: the
rows there labelled interior were whole-quiver rows, which made them stronger
claims about the *LNA* and weaker ones about *locality*, and the distinction
matters because locality is what a move rule is. `probe.py` reports the arrows a
run can rewrite for that reason, and lengthens the quiver unless told not to.

E-025.

---

## F-023 — What a rule needs to fire is a spectator, and the ends are where they are
*2026-09-16*

F-022 put sixteen rules against the ends of the quiver and coverage rose more
than the whole floating table had bought. Running discovery there properly says
how much more there is, and why the rules found until now so rarely fire.

**Discovery against an end.** `discover.py --anchor both --max-arrows 5
--max-width 6` plants each of the 74 patterns flush against each end of A_11 and
A_12, mutates within three vertices of it, keeps the rewrites described at both
lengths, and verifies each at the four lengths its window fits in:

| | |
|---|---|
| rewrites described | 1344, in 286 s |
| recurring at both lengths | 892 |
| verified, no failures | **724** |
| of those, a floating rule restricted to an end | 94 |
| genuinely anchored | **630**, 315 at each end |

315 at each end is a consistency check worth noticing: the relation dual
exchanges the two ends, so a rule at one has a mirror at the other, and the
counts had to come out equal.

**Most of them carry a spectator, and that is the finding.** A rule's window has
until now held nothing but the relations it rewrites -- `matchesAt` refuses a
position where any other relation reaches in. That is what makes a rule true and
it is why so few of them match anything: of the 155 LNAs left unplaced at n = 8,
**126 have a rule whose left-hand pattern is present and which does not fire**,
because one to four further relations are sitting in the window doing nothing.
Only 29 have no rule with their pattern at all. Among the 229 anchored rules
that change the orbit partition, **188 carry at least one relation that appears
unchanged on both sides** -- a bystander the rewrite steps around.

So the blockage was never that the patterns were too small. It was that they
were too clean.

**They do not compress into families.** Setting the spectators aside leaves 190
distinct rewrites among the 229, so unlike the slide families (F-013, F-020)
there is no statement covering many at once, and `quivermutation/endMoves.py`
lists them. 38 need one mutation, 99 two, 92 three.

**Two of them are worth reading on their own.**

```
window 2 arrows at the left end:   (0:2)  ->  -    via [1]
window 2 arrows at the right end:  (0:2)  ->  -    via [-3]
```

A lone relation of two arrows at an end of the quiver is simply deleted, by a
single mutation at the end vertex. That is operation 2 of
`cor:EquivNakayamaAlgebras` -- "a relation of two arrows does not change the
class" -- appearing as a mutation rather than as a theorem, and it is the rule
NOTES warned about: stated as a *floating* two-arrow window it holds 63 times
and fails 130, which is exactly right, because away from an end it is false.

**And widening the rules to admit a spectator is what the diagnosis was for.**
`lnaMoves.spectatorExtensions` puts one untouched relation into a rule's window,
growing the window by up to three arrows to make room, and `verifyMove` decides;
10609 such widenings produced **625** verified rules in eight minutes of no
searching at all, of which `spectatorMoves.SPECTATOR_MOVES` lists the 270 that
change the orbit partition -- 125 that float and 145 that need an end (E-024).

**What the two batches are worth.** LNAs placed with no mutation search:

| n | LNAs | theorem | + floating | + anchored | + widened |
|---|---|---|---|---|---|
| 6 | 42 | 34 (81%) | 35 (83%) | 42 (100%) | 42 (100%) |
| 7 | 132 | 89 (67%) | 95 (72%) | 127 (96%) | **132 (100%)** |
| 8 | 429 | 233 (54%) | 246 (57%) | 347 (81%) | **406 (95%)** |
| 9 | 1430 | 610 (43%) | 644 (45%) | 863 (60%) | **1038 (73%)** |

(E-025 has since raised the last column again, to 98% at n = 8 and 84% at n = 9,
and measured n = 10 and n = 11 for the first time: 63% and 47%.)

A classification of A_6 or A_7 is now a table lookup; A_8 needs a search for 23
rows and A_9 for 392. The rows still left are still exactly the heavily
overlapping ones, so F-021's reading is unchanged; there is simply much less of
it. And the diagnostic says the same thing about them as before: at n = 9, 359
of the 392 have a rule whose pattern is present and blocked by a bystander, and
only **33** have no rule with their pattern at all. Those 33 are what to look at
by hand.

**The curation, stated plainly.** `endMoves.DISCOVERED_END_MOVES` lists the 229
rules that change the orbit partition at n <= 9, not all 630 verified ones, and
`spectatorMoves.SPECTATOR_MOVES` the 270 of 625 the same way. A rule left out
reaches nothing the listed ones do not *at the lengths measured*, and could in
principle be the one that matters at n >= 10; re-running the commands gets them
all back. E-023, E-024.

---

## F-022 — An overlapping pair is frozen in the interior and comes apart at an end
*2026-09-16*

F-021 says the whole gap is the LNAs whose relations overlap in two or more
arrows, and that the commonest blocking configuration by a wide margin is a
**pair** of relations sharing two or more arrows. This is why the move rules
cannot place them, and what does.

**In the interior the overlap of an isolated pair does not move.** Plant
`(1:3) (2:3)` -- two relations of three arrows at consecutive vertices, sharing
two -- in the middle of A_13, with four arrows of empty quiver on the left and
six on the right, and enumerate every admissible mutation sequence near it:

| mutations | margin | LNAs reached | with a smaller maximum overlap |
|---|---|---|---|
| 3 | 3 | 8 | **0** |
| 4 | 3 | 14 | **0** |
| 5 | 3 | 22 | **0** |
| 4 | 6 | 34 | **0** |

The last row is stronger than the others and was first recorded as though it
were not: a margin of 6 around arrows 5 to 8 of A_13 admits *every vertex of the
quiver*, so it says that four mutations anywhere in A_13 -- ends included --
leave the pair intact.

Every LNA reachable still has two relations sharing two arrows. The same holds
for `(1:4) (2:4)` (overlap 3, 17 reached at three mutations, none below 3),
`(1:5) (2:5)` (overlap 4, 16 reached, none below 4) and for the unequal pairs
`(1:3) (2:4)` and `(1:4) (3:3)`, where the overlap goes *up* to 3 in three of
the sixteen but never down. Neither depth nor width is the obstacle: depth 5 and
a margin of 6 reach further into the quiver and find nothing new.

**A third heavily overlapping relation unlocks it; anything else does not.**

| pattern | run of | three mutations reach |
|---|---|---|
| `(1:3) (2:3)` | 2 | overlap 2, all 8 |
| `(1:3) (2:3) (3:3)` | 3 | **overlap 0**, via `[6, 5, 6]` |
| `(1:3) (2:4) (3:4)` | 3 | **overlap 0** |
| `(1:4) (2:4) (4:3)` | 3 | **overlap 0** |
| `(1:4) (2:4) (3:4)` | 3 | overlap 2, down from 3 |
| `(1:2) (2:3) (3:3)` | 2 | overlap 2, all 31 |
| `(1:3) (2:3) (4:2)` | 2 | overlap 2, all 31 |
| `(1:3) (2:3) (5:2)` | 2 | overlap 2, all 35 |

A third relation only helps when it *also* shares two or more arrows with the
pair: `(1:2) (2:3) (3:3)` has three relations and is as stuck as the bare pair,
because its first relation shares only one arrow. So the parameter is the length
of the **overlapping run** -- maximal relations linked by an overlap of two or
more, `overlap.overlapRuns` -- and a run of two is frozen where a run of three is
not.

**QUALIFIED 2026-09-16 by E-025.** "A run of three is not frozen" is true of the
short runs in the table above and false in general: `(1:3) (2:6) (3:7)` and
`(1:5) (2:6) (3:7)` reach two LNAs each at three mutations and neither lowers
the overlap. What the length of the run buys is probably mutations rather than
freedom -- F-020's one per arrow travelled -- so a long run may well dissolve
deeper down. Do not quote the sentence above without this. The rules that dissolve a run of three were already in the table; nothing in
it dissolves a run of two, and E-021 says why nothing was ever going to be found.

**At an end of the quiver the pair collapses in two mutations.** The source of
the line has no arrow into it, so a mutation there is not the mutation the same
rewrite would be in the interior. Where the pair starts at vertex 1, two *right*
mutations at vertex 1 delete the second relation outright; at the sink, two left
mutations at vertex n delete the first:

```
window l + 1 arrows at the left end:   (0:l) (1:l)  ->  (0:l)    via [1, 1]
window l + 1 arrows at the right end:  (0:l) (1:l)  ->  (1:l)    via [-(l+2), -(l+2)]
```

Verified for `l = 2` to 7, both ends, at the four lengths `l + 2 .. l + 5` each:
**9 confirmations apiece, no failures**, a confirmation being an admissible
sequence landing on the predicted LNA with the Coxeter polynomial kept. The
window is exactly the pair's span, so no other relation may touch it.
`lnaMoves.endPairCollapseRules` generates the family.

**Why this needed the framework to grow a notion it did not have.** Every rule
until now was a rewrite holding at *every* window position. The collapse holds at
one position and is false at all the others -- checked, not assumed:
`(0:3) (1:3) -> (0:3)` via `[1, 1]` stated as a floating rule fails, which is
`test_the_end_pair_collapse_is_false_in_the_interior`. So a description now
carries an optional anchor, `'left'` or `'right'`; `lnaMoves.windowStartsFor` is
the single gate every caller slides a rule through, and an anchored rule is
offered only its own position. Stating such a rule as though it floated is
exactly how R-009's false rules arose.

**What it buys, and it is more than sixteen rules should.** With the anchored
family and no other change:

| n | LNAs | theorem | + floating orbits | + anchored |
|---|---|---|---|---|
| 6 | 42 | 34 | 35 (83%) | **38 (90%)** |
| 7 | 132 | 89 | 95 (72%) | **107 (81%)** |
| 8 | 429 | 233 | 246 (57%) | **274 (64%)** |
| 9 | 1430 | 610 | 644 (45%) | **726 (51%)** |

Sixteen anchored rules place 82 rows at n = 9 where all 123 floating rules
placed 34. The mechanism is the pair slide (F-013) walking a pair to an end and
the collapse taking it from there -- `closureUnderMoves(10, [3,3,0,0,0,0,0,0])`
is the six positions of the pair *and* the two LNAs where it has lost a
relation, which is
`test_the_pair_slide_walks_a_pair_along_the_quiver_and_off_each_end`.

E-021, E-023. Tests: `tests/test_overlap.py`.

---

## F-021 — What a classification search still has to find is exactly the heavily overlapping LNAs
*2026-09-16*

H-003 asked what relation patterns the rows still needing a search actually
have. They have one, and it is sharp: **every LNA the quipu theorem and the move
orbits fail to place has two consecutive relations sharing two or more arrows,
and almost every LNA that has two such relations is one of them.**

**Overlap is the right coordinate because the theorem's condition is a bound on
it.** Consecutive relations `(n_i, l_i)`, `(n_{i+1}, l_{i+1})` share
`max(0, n_i + l_i - n_{i+1})` arrows, and *almost separate* -- the hypothesis of
`thm:QuipuToAn` -- is exactly that this never exceeds one. So "what the theorem
misses" and "what overlaps by two or more" are the same set, not merely
correlated ones; `test_almost_separate_is_exactly_overlap_at_most_one` checks the
two predicates against each other over every LNA of lengths 4 to 8.

**The measurement.** For each length, partition the LNAs into orbits under the
verified move rules and call an LNA *covered* when its orbit contains one the
theorem names -- which is what `seedTableFromQuipuTheorem` fills in before a
single mutation is computed. With the 123 floating rules:

| n | LNAs | overlap ≤ 1 | overlap ≥ 2 | of those, covered | left |
|---|---|---|---|---|---|
| 6 | 42 | 34 | 8 | 1 | 7 |
| 7 | 132 | 89 | 43 | 6 | 37 |
| 8 | 429 | 233 | 196 | 13 | 183 |
| 9 | 1430 | 610 | 820 | 34 | 786 |

Two things at once. Nothing at overlap ≤ 1 is ever left over -- so there is no
second phenomenon hiding among the rows the theorem does cover, and the search's
whole remaining job is the overlapping ones. And of the 820 heavily overlapping
LNAs at n = 9 the entire rule table reaches **34**: the rules found so far are
almost exactly the rules that keep an LNA where it already was.

**What the leftovers look like.** Their heavily overlapping runs -- maximal
relations linked by an overlap of two or more -- are overwhelmingly *pairs*, and
overwhelmingly the shortest pair there is:

| run | n = 7 | n = 8 | n = 9 |
|---|---|---|---|
| `(1:3) (2:3)` | 21 | 96 | 391 |
| `(1:4) (2:4)` | 6 | 40 | 198 |
| `(1:3) (2:4)` | 5 | 31 | 144 |
| `(1:4) (3:3)` | 5 | 31 | 144 |

and by the length of the longest run in the LNA, 434 of the 786 left at n = 9
have no run longer than two. That pair is what F-022 is about.

**The instrument.** `quivermutation/overlap.py` -- `overlapProfile`,
`maxOverlap`, `overlapRuns`, `coverage`, `blockingCores` -- and `overlaps.py`
over it:

```bash
python overlaps.py 6 7 8 9              # the table above
python overlaps.py 9 --cores            # what is left, by overlapping run
python overlaps.py 9 --floating         # without the rules anchored to an end
```

The whole thing is seconds per length because the orbits are computed as
rewrites on the relation lengths with no mutation run at all -- legitimate
exactly because every rule has been checked against the engine wherever it
applies (F-017).

H-003 → its question answered, its diagnosis refuted (R-010). E-020.

---

## F-020 — Three rule families whose mutation count grows with their parameter
*2026-09-15*

H-008 suspected that the move rules are members of families parameterised by
relation length and overlap, and that **a family can be simple to state while
needing more mutations for larger parameters** — which would make its later
members invisible to a search bounded at three. That is now a family, not a
suspicion.

**The lone short-relation slide.** A relation of two arrows with nothing else in
its window travels `d` arrows right under the `d` left mutations at the window's
vertices 3, 4, …, d + 2, and back under the `d` right mutations at d + 1, d, …, 2.
The window is `d + 2` arrows wide.

| d | window | sequence | found by |
|---|---|---|---|
| 1 | 3 | `[-3]` | E-010, two mutations |
| 2 | 4 | `[-3, -4]` | E-010 |
| 3 | 5 | `[-3, -4, -5]` | E-011, three mutations |
| 4 | 6 | `[-3, -4, -5, -6]` | **nothing** — needs four |
| 5 | 7 | `[-3, …, -7]` | needs five |
| 6 | 8 | `[-3, …, -8]` | needs six |
| 7 | 9 | `[-3, …, -9]` | needs seven |

Verified for d = 1 to 7, both directions, at the four lengths `d+3 .. d+6` each —
so up to A_14 and its 742900 LNAs: **63 confirmations apiece, no failures**, where
a confirmation is an admissible sequence landing on the predicted LNA with the
Coxeter polynomial kept. The count
being 63 at every d is itself a consistency check — the pattern is one relation
alone in its window, so the number of matching LNAs does not depend on how far it
travels.

**Two more of the same shape, read off consecutive widths in the enlarged table.**
Once the first family was recognised, looking for others cost minutes rather than
the hours a four-mutation search would have:

| family | rewrite on the window | sequence | mutations | verified |
|---|---|---|---|---|
| lone slide | `(0:2)` → `(d:2)` | `[-3, …, -(d+2)]` | d | d = 1–7, both directions |
| trailing walk | `(0:2) (2:2)` → `(1:2) (d+2:2)` | `[-3, -5, …, -(d+4)]` | d + 1 | d = 1–6 |
| spreading pair | `(1:2) (4:2)` → `(0:2) (d+4:2)` | `[2, -7, …, -(d+6)]` | d + 1 | d = 1–5 |

8 confirmations per member for the latter two, at three lengths each, no
failures; the spreading pair reaches A_14. Discovery had found d = 1 and 2 of each
and could not have found more — d = 3 of either needs four mutations. **One
mutation per arrow travelled** is the shape all three share: the relation that
moves furthest pays for each arrow, and any companion relation costs one more.

**Contrast with the pair slide (F-013).** That family is *two* mutations for every
relation length: only the window grows. Here the mutation count grows with the
parameter, so the two kinds are the two halves of H-008's statement, and this kind
explains why discovery keeps finding "new" rules that are the same rule.

**The practical consequence, and it is the point.** Discovery can only ever find
an **initial segment** of such a family — `d <= maxSteps` — so raising the search
bound by one buys one more member at multiplying cost, while recognising the
family gives every member at once. `shortRelationSlideRules`,
`trailingRelationWalkRules` and `spreadingPairRules` generate the three, as
`pairSlideRules` does for the other kind, taking the table from 96 listed rules
to 123. Before spending a four-mutation search, look at what the three-mutation
one found for a family whose members would be out of reach.

**Checked against the obvious objection, and it holds.** If a rule's mutation
count grew only because the search that found it looked at sequences running one
way, a sequence mixing left and right mutation might do the same job in fewer.
It does not. Of 150 rules of three mutations or more, three came back with a
shorter sequence -- all three members of this family, which is where a shortcut
would have mattered -- and none of the three is a rule: 1 confirmation against
21 failures as a floating rewrite, 1 against 8 anchored to an end. They hold at
the one LNA the search tried and nowhere else, for every d from 1 to 7 (E-026).

**What does not generalise.** The three slide families' inverses come for free:
reverse the sequence, negate each vertex, move each one step toward zero. That is
not a property of the table, though — it holds for only 12 of the 96 listed
rules, and fails or leaves the window for the rest (E-019). It works here because
a slide's sequence is one uniform run of mutations in one direction.

H-008 → CONFIRMED. E-011, E-019.

---

## F-019 — The gentle LNAs are one class, so gentle invariants separate nothing here
*2026-09-15*

Idea 22 (and NOTES idea 14) proposed the Avella-Alaminos–Geiss invariant as the
independent route to separating classes the quipu theorem does not reach, on the
stated ground that "LNAs are gentle". They are not, and even where they are the
route is empty.

**An LNA is gentle exactly when every relation has two arrows.** A gentle algebra
is a monomial algebra whose ideal is generated by paths of **length two**, plus
degree conditions on each vertex and two conditions on pairs of arrows. On a line
quiver every vertex has at most one arrow in and one out, so all of that is
automatic and only the length of the relations is at issue — and a relation of
three or more arrows is not generated by paths of length two. Such an LNA is
still a *string* algebra, but the AAG invariant is a gentle-algebra invariant.

**And every gentle LNA is the hereditary one.** Two facts collide:

1. all relations of length 2 implies almost separate relations — the condition
   `n_{i+1} >= n_i + l_i - 1` reads `n_{i+1} >= n_i + 1` when every `l_i` is 2,
   which distinct starts already give; and
2. operation 2 of `cor:EquivNakayamaAlgebras` drops a relation of two arrows
   without changing the derived equivalence class.

Drop them all and nothing is left. So the `2^(n-2)` gentle LNAs of length n are
all in the single class `P^(0)_(0,n-1)`, the class of the path algebra of A_n.
Checked against the classifications: 8, 16, 32, 64, 128 gentle LNAs at n = 5 to 9,
every one of them in that class and no other class containing any, and the
theorem-level statement pinned for n = 4 to 10 in
`test_every_gentle_lna_is_the_hereditary_one`.

**Consequence.** No gentle-algebra invariant can separate two LNA classes, because
no class but the hereditary one contains a gentle algebra at all. In particular
neither member of the cospectral pair `P^(1,4)_(1,0,1)` / `P^(1,2)_(1,1,2)` has a
single gentle member among its 18 — checked directly. The AAG route is closed;
what is left on idea 22's list is Hochschild cohomology, which is a derived
invariant of any algebra and does not care about the relation lengths. R-008.

---

## F-018 — A name read off the Coxeter polynomial must never separate two classes
*2026-09-15*

The n = 9 classification, re-run on the corrected engine (F-015) and the loosened
gate (F-016), came out at **22** classes rather than F-011's 20. The engine was
not at fault. The two extra classes were a **circular use of the Coxeter
polynomial** in the merge step, and once that is removed the partition is F-011's
20 again, on the very same table of rows.

**The mechanism.** A class the quipu theorem does not name gets a fallback name
from `canonicalWeightType`, which searches the weight types of the right order
for one whose canonical algebra has *exactly this class' Coxeter polynomial*. So
`C(2,3,5)` is a restatement of the polynomial and nothing more. But it is a
different *string* from `P^(5)_(1,2)`, and `mergeReport` separated two classes
whenever their form strings differed. The class named `C(2,3,5)` shared its
polynomial with `P^(5)_(1,2)` — which is how it got that name — and was declared
distinct from it on the strength of the very polynomial the two have in common.
The same happened to `C(2,2,6)` against `P^(1,1)_(1,3,1)`. Worse, the label also
stopped the class ever being searched again: `resolveMergeCandidates` only looked
at classes with *no* form.

**Both pairs really are one class, two independent ways.**

1. **Mutation.** `1;2;3;4;5|2;3;4;5;6;7`, of the class the run called `2334400`,
   reaches `1;2;3;4;5|2;3;4;5;6;7|6;7;8` of `P^(5)_(1,2)` in **two** mutations;
   `2;3;4;5|3;4;5;6|6;7;8;9`, of `2233030`, reaches `1;2;3;4|2;3;4;5|6;7;8;9` of
   `P^(1,1)_(1,3,1)` in **two**. Both also merge from the relation dual. These
   were available to the old resolution step at any depth it ran; it never looked.
2. **Theory.** A canonical algebra is derived equivalent to a hereditary algebra
   exactly when its weight type is *domestic* — `(p,q)`, `(2,2,n)`, `(2,3,3)`,
   `(2,3,4)`, `(2,3,5)` — with the partner of extended Dynkin type Ã, D̃, Ẽ6,
   Ẽ7, Ẽ8. Both (2,3,5) and (2,2,6) are domestic, and the trees come out at
   exactly the two quipus in question:

   | weight type | affine type | vertices | tree |
   |---|---|---|---|
   | (2,3,3) | Ẽ6 | 7 | `P^(2)_(2,2)` |
   | (2,3,4) | Ẽ7 | 8 | `P^(3)_(1,3)` |
   | (2,3,5) | Ẽ8 | 9 | `P^(5)_(1,2)` |
   | (2,2,6) | D̃8 | 9 | `P^(1,1)_(1,3,1)` |

   So a **domestic** `C(...)` can never be a class of its own: the class is a
   quipu class, and the quipu theorem has already named it under another name in
   the same table. A domestic weight type in that column is always a merge nobody
   found. `C(2,4,4)`, the n = 9 class that survives, is **tubular** — the boundary
   past which a canonical algebra is derived equivalent to no hereditary algebra
   — so that one is genuine.

**What the column means now.** Three kinds of value, and only two of them decide
anything:

* a quipu name or a tree encoding — **proved**, by a mutation path to a
  relation-free quiver or by the theorem. Two classes with different proved forms
  are different classes; two with the same one are the same class.
* `not piecewise hereditary` — **proved**, negatively. It separates such a class
  from every quipu class and merges nothing.
* `C(...)` — **not proved**. `isCoxeterDerivedForm` marks it; it may not separate
  and it may not merge, and a class carrying one stays a merge candidate.

**Why the pipeline also had to be reordered.** The cheap proof (a bounded search
for a path to an already-named class) now runs *before* the fallback that reads
the polynomial. Run the weak name first and it looks like an answer: a class one
mutation away from a quipu class gets labelled and is never searched again, which
is exactly what happened. `nameClassesFromTheorem` → `resolveMergeCandidates` →
`nameRemainingClasses`.

**One more thing this settles.** Every quipu of order n is realised by an LNA with
almost separate relations, so the quipu theorem names *every* quipu class in the
table. A class it leaves unnamed is therefore not a quipu class, and if it turns
out to be one after all, it must share a Coxeter polynomial with the named copy —
so it is always in a merge-candidate group, and the expensive search for a
relation-free quiver of its own (`--form-depth`) can never be what finds a quipu.
It is off by default for that reason.

E-017. Tests: `tests/test_merge_decisions.py`.

---

## F-017 — A move rule is a local rewrite, and the right encoding is per arrow
*2026-09-15*

H-009 suspects the move table is a one-dimensional cellular automaton, and names
its own caveat: a rule whose applicability depended on the row far away would not
be local at all, whatever it looked like. It does not. Both halves of that check
pass.

**Applicability is local.** `lnaMoves.matchesAt` scans the whole relation-length
row, but it does not need to. The same answer comes out of

  (a) the window's own cells -- a relation starting at an offset inside the
      window and staying inside it has at most `width` arrows, so what a cell
      inside can say is bounded by the width; and
  (b) **one bit**: whether any relation covers the window's first arrow having
      started strictly before it.

Checked exhaustively over every rule in `VERIFIED_MOVES`, every admissible LNA
and every window position at lengths 5 to 9: **991,064 comparisons, 1234 of them
an actual match, and no disagreement anywhere.**

**And the mutations a local match licenses are legal.** Re-verifying the whole
table where each rule fits gives **1218 confirmations and zero failures** --
every position where a rule matches, the mutation sequence is admissible at every
step, lands on the predicted relation lengths, and keeps the Coxeter polynomial.
So the admissibility side condition is not an extra non-local hypothesis riding
along; wherever the local pattern holds, the rewrite is legitimate.

**What the one bit says about the encoding, and this is the useful part.** Bit
(b) is exactly what the per-vertex encoding *cannot* supply locally: cell `i`
holds the number of arrows in the relation starting at vertex `i + 1`, and that
number is unbounded in `n`, so a relation can reach arbitrarily far to the right
and no fixed neighbourhood of cells sees it coming. Two consequences:

* the per-vertex row is **not** a good CA state -- unbounded alphabet, unbounded
  reach;
* a row indexed by **arrows** rather than vertices, each carrying whether it is
  covered and whether a relation starts or ends there, has a **fixed alphabet**
  and makes bit (b) a property of the cell at the window's edge. Over that
  encoding the move table is local with a margin of one cell.

So the CA reading is about the right object, provided the state is the arrow row.
That also says where the analogy will strain: converting a per-vertex row to a
per-arrow one needs to know how many relations are open at each arrow, which is
a counter, and the count is bounded only under *almost separate* relations
(overlap at most one arrow, so at most two). For the heavily overlapping LNAs --
which is exactly where the classification still needs a search (H-003) -- the
translation is not finite-state.

**Evidence.** `tests/test_lna_moves.py`,
`test_whether_a_move_applies_is_a_local_condition` at lengths 5 to 7 with 8 and
9 marked slow, and `test_each_rule_holds_wherever_it_applies` for the legality
half. E-016.

---

## F-016 — The search's gate is the paper's criterion now, and it costs nothing
*2026-09-15*

The admissibility gate was a **stricter** reading of the paper's theorem than the
theorem states, and it was applied twice:

* `mutation.mutationIsPossibleAtVertex` walked the relations and refused a vertex
  on any minimal zero relation whose last arrow left it and whose truncation was
  not itself written as a relation;
* `search.mutationSearchDepthFirst` then counted, for every predecessor `v` and
  every arrow `i -> w`, the paths `v -> i` against the paths `v -> w`, and
  refused the vertex if any one arrow lost a path.

Both are the same idea and both are "**every** arrow out of the vertex must keep
every nonzero path nonzero". The paper's theorem rules mutation out only when a
nonzero path ending at the vertex dies against **every** arrow out of it — so
where the two differ, the old gate was refusing mutations the paper allows. Both
also decided "nonzero" syntactically: a zero relation written inside the path, or
a path count taken up to the commutativity relations.

The gate is now `procedure.isMutable`, the theorem read exactly with "nonzero"
decided over the ideal, and it is the only gate — the inline path count is gone.

**What it changed.**

| | |
|---|---|
| mutations newly allowed, n = 5 and 6 to depth 3, n = 7 to depth 2 | **280** |
| of those, Coxeter polynomial preserved | **280** — all of them |
| mutations the new gate refuses that the old one allowed | **0** |
| n = 6 classification | same 4 classes, same sizes |
| n = 7 classification | same 6 classes, same sizes |
| n = 8 classification | same 11 classes, same sizes |
| n = 7 classification time | **38 s → 20 s** |

Faster, which is the opposite of what a more permissive gate suggests: the old
gate cost two syntactic sweeps per vertex — `allRelsInPathAlgebra` at every
search node, then a path count per (predecessor, successor) pair — where the new
one answers from the ideal directly. Removing the sweeps more than pays for the
extra branches.

**Why every newly allowed mutation had to be checked.** The criterion is
*necessary and not sufficient*: the theorem's hypothesis is
`Hom(P_i*[1], Lambda) = 0`, and the paper says explicitly that this is in
general not equivalent to a condition on the quiver. So a mutation the gate
newly allows could in principle fail to be a derived equivalence, and the
Coxeter polynomial would move across it. None does. This is the check R-005
exists to insist on, and being stricter than the paper was the old gate's way of
avoiding having to make it.

**What did not change.** `A_{7,(2,4)}^{(3,3)}` still reaches no relation-free
quiver within depth 8 — the loosened gate does not rescue it, so the point of
`test_the_theorem_answers_where_the_search_gives_up` stands: there are classes
the theorem names outright that no affordable search reaches.

**Evidence.** `tests/test_procedure.py` holds the old criterion, as
`strictlyMutable`, and checks the two against each other: they agree on every
vertex of every LNA of lengths 4 to 6 (a line's vertices have one arrow out,
which is where they coincide), and after one mutation every disagreement is in
the permissive direction, at a branching vertex, with the Coxeter polynomial
preserved. The wider sweep is E-015.

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

*Amended 2026-09-15: confirmed on the corrected engine, after F-018. The first
re-run gave 22, from a circular use of the Coxeter polynomial in the merge step
and not from the engine; the two extra classes were `C(2,3,5)` and `C(2,2,6)`,
which are the domestic weight types of `P^(5)_(1,2)` and `P^(1,1)_(1,3,1)`. The
count and the partition below are unchanged.*

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
