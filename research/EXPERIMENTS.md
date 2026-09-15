# Experiments

Runs made, newest first, with parameters and outcome — including runs that found
nothing, which are recorded precisely so they are not repeated. See
[`README.md`](README.md).

---

## E-018 — How far the gentle condition reaches
*2026-09-15* · **one class, the hereditary one** → F-019, R-008

Idea 22's premise, checked before implementing anything.

Every LNA of lengths 5 to 9 tested for `all(arrows in (0, 2))`, then grouped by
the class the classification puts it in:

| n | gentle LNAs | 2^(n-2) | classes containing one |
|---|---|---|---|
| 5 | 8 | 8 | `P^(0)_(0,4)` only |
| 6 | 16 | 16 | `P^(0)_(0,5)` only |
| 7 | 32 | 32 | `P^(0)_(0,6)` only |
| 8 | 64 | 64 | `P^(0)_(0,7)` only |
| 9 | 128 | 128 | `P^(0)_(0,8)` only |

So every gentle LNA is in the class of the path algebra of A_n, and the other 21
classes at n = 9 contain none -- including both members of the cospectral pair,
which have 18 members each and not one gentle among them. Seconds to run, against
the days an AAG implementation would have taken.

The reason is a theorem, not a coincidence: all relations of length 2 implies
almost separate relations, and operation 2 of `cor:EquivNakayamaAlgebras` drops
such a relation without changing the class, so dropping them all leaves the
hereditary algebra. Pinned for n = 4 to 10 in
`test_every_gentle_lna_is_the_hereditary_one`.

Also established while looking: `WebSearch` reaches the literature from the
session sandbox even though `curl` and `WebFetch` to arxiv.org are refused by the
egress proxy. Enough to find and identify a paper, not to read one.

---

## E-017 — n = 9 on the corrected engine
*2026-09-15* · **22 classes, then 20** → F-018

`classify.py 9` on the engine of F-015 with the gate of F-016, default depths
(`--depth 6 --resolve-depth 6`). About 90 minutes.

The run placed all 1430 rows and left nothing a candidate, but at **22** classes:
18 quipus plus `C(2,3,5)` (46 LNAs), `C(2,2,6)` (13), `C(2,4,4)` (8) and one not
piecewise hereditary. It reported three groups "proved distinct despite sharing a
Coxeter polynomial", two of which were `C(2,3,5)` against `P^(5)_(1,2)` and
`C(2,2,6)` against `P^(1,1)_(1,3,1)`.

**Diagnosis, in this order.**

1. The `C(...)` names come from `canonicalWeightType`, which reads them off the
   class' Coxeter polynomial — so they cannot separate two classes that share
   one. Circular.
2. A direct probe: iterative deepening from every member of each of the two
   classes, and from every member's relation dual, reporting every foreign class
   reached. Both merged **at depth 2**, from the first member tried and from its
   dual as well. Seconds, against the 90 minutes of the run.
3. (2,3,5) and (2,2,6) are domestic weight types, and their extended Dynkin trees
   are `P^(5)_(1,2)` and `P^(1,1)_(1,3,1)` — the two classes they were separated
   from. Computed, not asserted.

**Re-run of the post-search half only** (the rows were sound; only the merge step
was wrong), on the same table: both merged at the first depth tried, then
`3033030` certified not piecewise hereditary and `3345000` named tubular
`C(2,4,4)`. **1430 LNAs, 20 classes, 0 candidates, 1 separated** — the separated
group being the cospectral pair of F-010, exactly F-011. Under a minute.

A whole run from scratch on the fixed pipeline is in flight; this entry gets
its result when it lands.

---

## E-016 — Are the move rules local?
*2026-09-15* · **yes, both halves** → F-017

H-009's own caveat, checked before anything else was built on it.

1. **Applicability.** `matchesAt` against a predicate reading only the window's
   cells plus one bit (a relation covering the window's first arrow having
   started earlier), over every rule in `VERIFIED_MOVES` x every admissible LNA
   x every window position:

   | n | comparisons | matches | disagreements |
   |---|---|---|---|
   | 5, 6, 7 | 67,712 | 150 | 0 |
   | 8, 9 | 924,352 | 1084 | 0 |

2. **Legality.** The whole table re-verified where each rule fits: **1218
   confirmations, zero failures** -- every match is an admissible sequence
   landing on the predicted LNA with the Coxeter polynomial kept.

Cheap: seconds for lengths 5 to 7, a couple of minutes for 8 and 9, and about
four minutes for the legality half. Rows 1 and 2 are tests now
(`test_whether_a_move_applies_is_a_local_condition`,
`test_each_rule_holds_wherever_it_applies`), so **do not repeat them by hand**.

The result that was not the question: the state has to be the **arrow** row, not
the vertex row -- see F-017. Anyone starting the CA literature sweep should start
there rather than from `relLengths`.

---

## E-015 — Every mutation the loosened gate newly allows
*2026-09-15* · **280 of them, all Coxeter-preserving** → F-016

Before switching the search's gate from the strict reading to the paper's
criterion, every mutation the switch would newly allow was enumerated and
checked. Walking out of every LNA of the length, at every vertex of every quiver
reached, comparing the old gate against the new one and computing the Coxeter
polynomial wherever they disagreed:

| n | depth | allowed by both | newly allowed | Coxeter moved | new gate narrower |
|---|---|---|---|---|---|
| 5 | 3 | 304 | 22 | 0 | 0 |
| 6 | 3 | 1450 | 138 | 0 | 0 |
| 7 | 2 | 1938 | 120 | 0 | 0 |

The last column matters as much as the others: a criterion that was *narrower*
anywhere would have meant the switch loses a mutation the published runs used,
and it never is.

Then the classifications, which are the acceptance test: n = 6, 7 and 8 all give
the same classes with the same sizes as before, and n = 7 dropped from 38
seconds to 20.

The old criterion is kept in `tests/test_procedure.py` as `strictlyMutable`,
which is what makes the comparison re-runnable; the first two rows are a test
now. **Do not repeat the n = 7 row** — about four minutes, and it says the same
thing as the other two.

---

## E-014 — The procedure on coefficients, against the one it replaced
*2026-09-14* · **agreement everywhere but two cases, which are R-007** → F-015

Five runs, all gated on `mutationIsPossibleAtVertex` so both implementations walk
the same mutations:

1. **One mutation, every admissible vertex, every LNA of n = 4..8.** 45 + 126 +
   462 + 1716 = 2349 comparisons, **zero** differences. Minutes.
2. **Depth-3 walks, n = 5 and 6.** 1446 and 7496 step comparisons, **zero**
   differences.
3. **Depth-3 walks, n = 7.** 37470 step comparisons, **2** differences, both
   after three mutations, both a relation the old implementation did not
   produce. These are the whole of R-007.
4. **The exact cleanup on the old steps' output**, n = 5 and 6 at depth 3 and
   n = 7 at depth 2: 1446 + 7496 + 4710 = 13652 comparisons, **zero**
   differences. Worth having separately, because it says the disagreement is in
   step 7 and not in the cleanup.
5. **Coefficients against the guess**, over every quiver within depth 3 of every
   LNA of n = 5 and 6: 1239 Cartan matrices, **zero** differences.

**A mixed engine is not an option, and this is how that was learned.** Running
the old steps 1-7 with the exact cleanup passed run 4 above and then reached
*two* different hereditary forms from `A_6` `3030` at depth 7 — a degree-4 tree
alongside `P^(1,1)_(1,0,1)`, which cannot both be one class. The exact cleanup
expects the relations step 7 produces; with step 7's output missing a relation it
cuts the wrong generators. Use one engine or the other, whole.

**Timings**, n = 7 over the 462 admissible single mutations: procedure 0.26 s
against 0.94 s, admissibility 0.14 s against 1.74 s. The exact versions are
3.6x and 12x *faster*.

**Do not repeat runs 1, 2, 4 and 5** — they are `tests/test_procedure.py` now.
Run 3 at n = 7 depth 3 takes about eight minutes and is worth re-running only if
step 7 changes.

---

## E-013 — Audit of the quipu symmetry, after R-006 was challenged
*2026-09-14* · **no defect found** → F-014

Four runs, in increasing cost:

1. **Canonicalisation against `networkx.is_isomorphic`**, over every quipu
   parameter pair of orders 3–11 (12 names at order 3 up to 28656 at order 11).
   Same canonical parameters iff isomorphic graphs, both directions, zero
   exceptions. Seconds.
2. **The paper's class-preserving operations against the quipu fibres**, over
   every LNA of lengths 4–10 with almost separate relations and no length-2
   relation. Orbits equal fibres exactly at every length; largest orbit 8, the
   paper's bound. Seconds. Now `tests/test_quipu_symmetry.py`.
3. **Tree enumeration against parameter enumeration.** The old
   `generateAllQuipus` (enumerate non-isomorphic trees, test the degrees) and
   `quipuForms.allQuipusOfOrder` (enumerate the P^(m)_(k) parameters,
   canonicalise) agree on the counts for orders 4–12, and no tree the first
   accepts is rejected by `quipuForms.isQuipu`. The old function only tested the
   "degree-3 vertices lie on one path" condition when there were more than three
   of them, which looked like a hole — but three branch vertices in a tree of
   maximum degree 3 always do lie on one path, since a path through two of them
   passes through the third, so four is the smallest number that can fail.

   Kept, since it is a genuinely independent route: it is now
   `quipuForms.quipusByTreeEnumeration`, with the degree test written out as
   `isQuipuByDegrees`, and the agreement is a test rather than a note here.
4. **Hereditary form by mutation search**, from all four long-relation members of
   `P^(1,4)_(1,0,1)` (`0003030`, `3030000`, `3060000`, `6000030`) and both of
   `P^(1,2)_(1,1,2)` (`0400030`, `3004000`), at depth 6, plus iterative
   deepening 2–6 from `3060000` and `3004000`. **Nothing reached** — no
   relation-free quiver from any of them. Tens of minutes.

Run 4 is the one that would have been independent of `thm:QuipuToAn`, and it is
simply out of range here, the same way `A_{7,(2,4)}^{(3,3)}` is (see the test
`test_the_theorem_answers_where_the_search_gives_up`). **Do not repeat it at
depth 6 or less.** Depth 7+ at n = 9 was not attempted and is expected to be
hours; the cheaper route to an independent check is a derived invariant computed
from the algebra, not a deeper search.

---

## E-012 — Pair slide at relation lengths 2 to 7
*2026-09-14* · **confirmed a family**

`lnaMoves.verifyMove` on the pair-slide rewrite for each `l`, both directions,
over lengths `l+3 .. l+6`. 22 confirmations per direction per length, zero
failures throughout. → F-013, confirming H-001.

Seconds to run. Should have been the first thing tried after finding the rule at
`l = 3`.

---

## E-011 — Interior discovery, three mutations
*2026-09-14* · **running**

`lnaMoves.discoverLocalMoves`, 26 patterns of up to 3 relations spanning ≤ 5
arrows, planted at offset 4 in A_13 and offset 5 in A_14, `maxSteps=3`,
`margin=3`, then verified over lengths 7..10.

    python interior.py 3 3 5

Cost: roughly 30 s per (pattern, embedding) at `maxSteps=3`. Tests H-007.

---

## E-010 — Whole-quiver discovery at length 8
*2026-09-13* · **34 rules**

`discoverMoves([8], maxSteps=2)`: 101 candidates, 85 not already known, 34
verified. All of window 6 — the width that first has room to sit clear of both
ends at that length.

Each was confirmed only 3 times at lengths 7–8, which is thin, so all 34 were
**re-verified over lengths 7 to 10**: all survived, 22 confirmations each. Keep
doing this for wide rules found at short lengths.

---

## E-009 — Whole-quiver discovery at length 7, three mutations
*2026-09-13* · **2 rules** — poor yield

`discoverMoves([7], maxSteps=3)`: 61 candidates, 45 new, **2** verified.

The yield is low because `describeLink` only admits a *local* rewrite — one whose
window contains every relation it touches — and at length 7 a three-mutation
sequence usually disturbs the whole quiver, so nothing recurs across positions.
This is the experiment that motivated interior embedding (H-007). **Do not repeat
at this length.**

---

## E-008 — Classification of n = 10
*2026-09-13, updated 2026-09-14* · **unfinished — resume it**

`classifyLength(10)`, several attempts, none yet complete. Furthest reached:
about 1900 of 4862 rows.

**Long runs do not survive.** Three separate causes, all worth knowing:

1. two attempts were killed by over-broad `pkill -f` patterns issued by the
   session itself — a pattern that also matches the shell issuing it kills the
   shell, and anything sharing its process group;
2. one died under `setsid` when the machine went away between sittings;
3. n = 10 takes hours, so any of the above is likely to happen at least once.

**So: resume rather than restart.** The table is written after every class
searched, and `classify.py --resume` continues from the existing CSV:

    python classify.py 10 --resume

Before assuming a long job is still running, check `ps` — a stalled row count
looks the same as a dead process.

---

## E-007 — Certificate propagation by vertex deletion
*2026-09-13* · **0 / 0 / 1 / 24 / 308**

`notPiecewiseHereditaryByDeletion` over every LNA of lengths 4 to 11 → F-012.
The zero below length 9 is the correctness check, not an absence of result.

---

## E-006 — Cospectral quipu enumeration to order 13
*2026-09-13* · **the collision map**

`quipuForms.cospectralQuipuGroups(n)` for n = 4..13, cross-checked against equal
Coxeter polynomials computed through each algebra's Cartan matrix for n = 4..11.
The two agree exactly. → F-010.

Seconds to run, no mutation search involved. `python classify.py <n> --collisions`.

---

## E-005 — Orbit verification under the move table
*2026-09-13* · **clean**

Every LNA of lengths 5 to 9: compute its orbit under the verified moves, apply
each recorded mutation sequence to check it reaches the class it claims, and check
the Coxeter polynomial is constant on the orbit. 1764 orbit members, zero
failures.

Run **before** trusting any change to the rule table — it is what caught R-005.

---

## E-004 — Reduction preserves the Cartan matrix
*2026-09-13* · **clean, after R-003**

Every legal mutation of depth ≤ 3 out of every LNA of lengths 5 to 8: 38095
reductions, zero changes. Exact and heuristic Cartan matrices agree throughout.
→ F-008.

First run gave 9 apparent failures; all were R-003, not the reduction.

---

## E-003 — Exact against heuristic Cartan matrix
*2026-09-13* · **agree everywhere tested**

All 624 LNAs of length ≤ 8, and all 8101 quivers reached by walking every legal
mutation of depth ≤ 3 out of all 188 LNAs of lengths 5 to 7. No disagreement, so
no published Coxeter polynomial moves. The shapes where the two models differ
(F-004) have not turned up in an LNA search.

---

## E-002 — Classifications of n = 5 to 9
*2026-09-12 – 2026-09-13* · **match the published table**

| n | LNAs | classes | time |
|---|---|---|---|
| 6 | 42 | 4 | ~13 s |
| 7 | 132 | 6 | ~37 s |
| 8 | 429 | 11 | ~4 min |
| 9 | 1430 | 20 | ~56 min |

n = 6, 7, 8 match arXiv:2305.06642 exactly. n = 9 → F-011. Lengths 6–8 are pinned
as `slow` tests.

---

## E-001 — Reproducing the published n ≤ 8 classification
*2026-09-12* · **the baseline**

The first cross-check of the restored code against the papers: relation-set counts
against the Catalan numbers, the worked example of arXiv:2112.08129 step by step,
Coxeter polynomials of A_n and D_n, and the class membership of the n ≤ 8 table.
Everything agreed once F-001 was fixed.
