# Experiments

Runs made, newest first, with parameters and outcome — including runs that found
nothing, which are recorded precisely so they are not repeated. See
[`README.md`](README.md).

---

## E-024 — Widening the rules to tolerate a bystander
*2026-09-16* · **625 verified, and A_7 needs no search at all** → F-023, H-011

F-023 said what stops a known rule from firing is a relation in its window that
it does not touch. This is that read as a construction rather than a diagnosis:
take each rule in the table, put one untouched relation -- a *spectator* --
somewhere in its window, growing the window by up to three arrows to make room,
and let `verifyMove` decide.

```bash
python discover.py --extend --jobs 4
```

| stage | |
|---|---|
| widenings generated from the 368 rules then in the table | 10609, in 13 s |
| of those, firing on an LNA no search had placed at n <= 9 | **1008** |
| verified with no failures | **625**, in 474 s |
| changing the orbit partition at n <= 9 | **270** -- 125 floating, 145 anchored |

Eight minutes, no search, no mutation run speculatively. The filter is what makes
it affordable and is worth keeping: generate freely, throw away everything that
would not fire on a row still needing a search, then verify.

**Coverage with no search at all.**

| n | theorem | before this batch | after |
|---|---|---|---|
| 6 | 81% | 100% | 100% |
| 7 | 67% | 96% | **100%** |
| 8 | 54% | 81% | **95%** -- 23 rows left of 429 |
| 9 | 43% | 60% | **73%** -- 392 of 1430 |

A classification of A_7 is now a table lookup. At n = 8 twenty-three rows need a
mutation search and at n = 9, 392.

**The number H-011 said to watch stayed small.** Re-running the blocked-rule
diagnostic -- for each unplaced LNA, the rule whose left-hand pattern is present
with the fewest extra relations in its window:

| | n = 8 before | n = 8 after | n = 9 after |
|---|---|---|---|
| blocked by a bystander | 126 | 18 | 359 |
| **no rule has this pattern at all** | 29 | **5** | **33** |

So the mechanical half is still the whole story: what is left is overwhelmingly
more of the same, and another widening pass is the obvious next run.

**And the 33 turn out to say the same thing.** Looked at by hand, they are
almost all a pair of relations one of which is *long*: `(1:3) (2:6)`,
`(1:5) (2:6)`, `(1:5) (2:7)` and the like. 32 of the 33 contain a relation of
five arrows or more and 23 contain one of six or more -- and discovery has never
been given a pattern like that. Both runs used `--max-arrows 5 --max-width 6`,
so a six-arrow relation could not appear beside another at all. This is F-013's
lesson again: *absence of a pattern from the table is evidence about the search,
not about the mathematics*. Raise the bounds before concluding anything about
these.

**What this does not say.** Every one of these rules was verified, but 355 of
the 625 are not listed, and the 270 that are were chosen for changing the
partition at n <= 9. That is a curation against the lengths measured, not a
claim about n >= 10. Re-run the command.

---

## E-023 — Discovery against the ends of the quiver
*2026-09-16* · **630 anchored rules, and coverage at n = 6 becomes complete** → F-023

The first run of `discoverAnchoredMoves`, on the framework F-022 added.

```bash
python discover.py --anchor both --max-arrows 5 --max-width 6 --jobs 4 --verify-cap 12
```

74 patterns x 2 ends x 2 lengths (A_11 and A_12) = 296 searches at three
mutations, margin 3.

| stage | |
|---|---|
| rewrites described | 1344, in 286 s |
| recurring at both lengths | 892 |
| verified with no failures | 724, in 70 s |
| a floating rule restricted to an end | 94 |
| genuinely anchored | **630** -- 315 at each end |
| changing the orbit partition at n <= 9 | **229**, and those are what is listed |

**A first attempt at the same run had to be abandoned**, and why is worth
recording. `verifyMove` enumerated the LNAs by building a path algebra for each,
which at length 12 is 58786 of them and six seconds -- per rule, and there were
886 to check. The enumeration is now cached as relation-length rows
(`nakayama.allRelationLengths`), the algebras built only where a rule actually
matches: 0.3 s instead of 6.5, and the verification of all 886 fell from hours
to 70 seconds. Anything that verifies many rules over the same lengths should go
through that function.

**What the run cost and bought.** Ten minutes end to end. Coverage with no
search at all: 100% at n = 6 (from 83%), 96% at n = 7 (72%), 81% at n = 8 (57%),
60% at n = 9 (45%).

**What is still not reached, and the next question.** 567 rows at n = 9, all of
them heavily overlapping, 306 at overlap 2. The diagnostic that pointed at the
spectators -- for each unplaced LNA, the rule whose left-hand pattern is present
with the fewest extra relations in the window -- says at n = 8 that 126 of 155
are blocked by a bystander and 29 by having no rule at all. Run it again after
the next batch: the count of "no rule has this pattern" is the one to watch,
because it is the part more discovery cannot fix.

---

## E-022 — Whether the relation dual widens the move orbits
*2026-09-16* · **it halves the orbit count and adds no coverage**

The relation dual -- reverse every arrow, renumber -- is one of the three
class-preserving operations of arXiv:2305.06642 and holds for *any* LNA, not
only an almost separate one (`nakayama.relationDual`). It is free, it is not in
the move orbit, and the obvious thought is that adding it would carry rows
across the overlap line for nothing. It does not.

| n | floating | + anchored | + anchored + dual |
|---|---|---|---|
| 7 | 95 covered, 84 orbits | 107, 71 | 107, **45** |
| 8 | 246, 310 | 274, 277 | 274, **156** |
| 9 | 644, 1106 | 726, 1019 | 726, **542** |

The orbit count roughly halves at every length and the covered count does not
move by one row. The reason is structural rather than accidental: the almost
separate condition is itself dual-symmetric, so the dual maps seeded to seeded,
and the rule table already contains the mirror of every rule it contains, so it
maps orbit to orbit. The dual therefore identifies orbits pairwise and never
joins a covered one to an uncovered one.

**Worth knowing, and worth not repeating.** Halving the orbit count is real and
would be worth having if orbits were the expensive object; they are not, the
uncovered rows are. Do not reach for the dual again expecting coverage.

---

## E-021 — What can be done to a heavily overlapping run, in the interior
*2026-09-16* · **the pair is frozen, a run of three is not** → F-022, R-010

The experiment H-003 asked for, aimed where F-021 says to aim it. Each pattern
planted in the middle of A_13 at offset 4 -- four arrows of empty quiver on the
left, six on the right -- with `lnaMoves.localMutationSequences` enumerating
every admissible sequence at vertices within the margin, and the reached LNAs
reported by maximum overlap.

**The isolated pair, at four settings.**

| pattern | start overlap | mutations | margin | reached | any lower |
|---|---|---|---|---|---|
| `(1:3) (2:3)` | 2 | 3 | 3 | 8 | no |
| `(1:3) (2:3)` | 2 | 4 | 3 | 14 | no |
| `(1:3) (2:3)` | 2 | 5 | 3 | 22 | no |
| `(1:3) (2:3)` | 2 | 4 | 6 | 34 | no |
| `(1:4) (2:4)` | 3 | 3 | 3 | 17 | no |
| `(1:5) (2:5)` | 4 | 3 | 3 | 16 | no |
| `(1:3) (2:4)` | 2 | 3 | 3 | 16 | no (3 of them go **up** to 3) |
| `(1:4) (3:3)` | 2 | 3 | 3 | 16 | no (3 go up to 3) |

About a quarter of an hour in total, the depth-5 probe a third of it. Neither depth nor margin
is the dial: doubling the margin at four mutations reaches 34 LNAs instead of
14 and not one of them has a smaller overlap.

**The margin-6 row is stronger than it was written as, and the description was
wrong.** A margin of 6 around a pattern at arrows 5 to 8 of A_13 admits the
vertices 1 to 13 -- *every vertex of the quiver*, both ends included. So that
row is not a probe of the interior at all: it says that from `00003300000`,
**four mutations anywhere in A_13** reach 34 LNAs and none of them has a smaller
overlap. That is a claim about the LNA rather than about locality, and it is the
stronger one. It was recorded here as an interior probe with a wide margin,
which it was not. `probe.py --allow-ends` is how to ask that question on
purpose; without the flag the quiver is lengthened to keep the ends out of
reach, so an interior probe stays one.

**A third relation, and which third relations count.**

| pattern | overlapping run | three mutations |
|---|---|---|
| `(1:3) (2:3) (3:3)` | 3 | down to **0**, via `[6, 5, 6]` |
| `(1:3) (2:4) (3:4)` | 3 | down to **0** |
| `(1:4) (2:4) (4:3)` | 3 | down to **0** |
| `(1:4) (2:4) (3:4)` | 3 | down to 2, from 3 |
| `(1:2) (2:3) (3:3)` | 2 | 31 reached, none lower |
| `(1:3) (2:3) (4:2)` | 2 | 31 reached, none lower |
| `(1:3) (2:3) (5:2)` | 2 | 35 reached, none lower |

The parameter is the length of the run of relations linked by an overlap of two
or more, not the number of relations present: `(1:2) (2:3) (3:3)` has three
relations and is as frozen as the bare pair, because its first shares one arrow
and not two. The rewrites that dissolve a run of three were already in the
table, found at length 8 and in E-011 -- so nothing here is a new rule, and that
is the result. **Do not run a deeper interior search for a rule that pulls an
isolated pair apart**; four probes at three settings of depth and two of margin
say there is none to find, and F-022 says where the pair does come apart.

Reproduce with `lnaMoves.localMutationSequences(13, relLengths, lo, hi, steps,
margin)` on `lnaMoves.embedPattern(13, pattern, 4)`.

---

## E-020 — What the theorem and the move orbits reach, by relation overlap
*2026-09-16* · **the gap is exactly overlap two and above** → F-021

The measurement H-003 has been asking for since 2026-09-13, now that there is a
coordinate to make it in. Every LNA of a length partitioned into orbits under
the verified rules -- applied as rewrites on the relation lengths, with no
mutation computed, which F-017 licenses -- and an orbit called covered when it
contains one the quipu theorem names.

With the 123 floating rules:

| n | LNAs | overlap 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|---|
| 6 | 42 | 16/16 | 18/18 | 1/7 | 0/1 | | | |
| 7 | 132 | 32/32 | 57/57 | 4/33 | 2/9 | 0/1 | | |
| 8 | 429 | 64/64 | 169/169 | 9/132 | 4/52 | 0/11 | 0/1 | |
| 9 | 1430 | 128/128 | 482/482 | 24/484 | 10/247 | 0/75 | 0/13 | 0/1 |

covered over total at each maximum overlap. Two readings, and both matter.
Everything at overlap 0 or 1 is covered, at every length -- which is the almost
separate set exactly, so the theorem's reach is not merely *mostly* the low
overlap rows, it is precisely them. And above the line the table reaches 34 rows
out of 820 at n = 9, none at all past overlap 3.

The leftovers' heavily overlapping runs, commonest first at n = 9: `(1:3) (2:3)`
391 times, `(1:4) (2:4)` 198, `(1:3) (2:4)` and `(1:4) (3:3)` 144 each. By the
longest run in the LNA, 434 of the 786 have nothing longer than a pair.

Seconds per length. `python overlaps.py 6 7 8 9 --cores`, and the same numbers
are pinned in `tests/test_overlap.py`.

---

## E-019 — Two more families, and whether a rule's inverse is free
*2026-09-15* · **two families confirmed, the inverse shortcut refuted** → F-020

Acting on F-020's own lesson rather than raising the search bound.

**The families.** Two candidates read straight off consecutive window widths in
the enlarged table, then verified at three lengths each with `verifyMove`:

| family | d | result |
|---|---|---|
| `(0:2) (2:2)` → `(1:2) (d+2:2)` via `[-3, -5, …]` | 1–6 | 8 confirmations each, no failures |
| `(1:2) (4:2)` → `(0:2) (d+4:2)` via `[2, -7, …]` | 1–5 | 8 confirmations each, no failures |

Discovery had found d = 1 and 2 of each; d = 3 needs four mutations and d = 6
needs seven, so the rest were out of reach of any search run so far. Minutes to
check, against the hours a four-mutation run costs. Both are generated now.

**The inverse shortcut, and it does not hold.** Family A's listed left slide is
exactly its right slide with the sequence reversed, each vertex negated, and each
then moved one step toward zero — which looked like it might be a property of the
window's numbering and so give every rule's inverse for nothing. Applied to all
96 listed rules and verified:

| sequence | inverts | fails | inverse leaves the window |
|---|---|---|---|
| all one direction | 12 | 20 | 28 |
| mixed directions | 0 | 5 | 31 |

So **12 of 96**. The transform works for the slide families because a slide's
sequence is a single uniform run; it is not a general fact about the table, and
the spreading pair -- whose sequence mixes a right mutation with left ones -- is
the counterexample closest to hand. Do not try it again.

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

**Then re-run whole, from scratch, on the fixed pipeline: the same 20**, with the
same sizes class for class. The naming order is visible in the log -- the theorem
names 18 classes, `resolveMergeCandidates` then merges 11 away (`2233030` into
`P^(1,1)_(1,3,1)` and `2334400` into `P^(5)_(1,2)` among them, both at depth 6),
and only then do the fallbacks name what is left: `3033030` not piecewise
hereditary by Proposition A9, `3345000` the tubular `C(2,4,4)`. One separated
group, the cospectral pair of F-010. About 35 minutes, 10 of them the search.

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
*2026-09-14, concluded 2026-09-15* · **44 rules, and 30 false ones caught** → F-020, R-009

`lnaMoves.discoverLocalMoves`, 26 patterns of up to 3 relations spanning ≤ 5
arrows, planted at offset 4 in A_13 and offset 5 in A_14, `maxSteps=3`,
`margin=3`. Re-run as

    python discover.py --jobs 2

after the original `interior.py` turned out never to have been committed.

**Discovery.** 52 searches, 315 s on two cores. 336 rewrites described, **166
recurring across both embeddings**, 134 of them not already in the table.

**Verification, first attempt — wrong, and instructively so.** All 134 checked at
the fixed lengths 7 to 10, which E-010 had used: 74 passed. But the lengths have
to follow the window, and a window of 9 arrows fits in A_10 at exactly one
position, flush against both ends. Each window-9 rule therefore got one
confirmation from one length.

**Verification, redone per rule at `width + 1 .. width + 4`.** 44 survive; **all
30 window-9 rules fail at length 11**, where the window can sit clear of the
ends — wrong rules, not thin ones. R-009.

| window | rules | lengths checked | confirmations |
|---|---|---|---|
| 5 | 2 | 6, 7, 8, 9 | 22 |
| 6 | 16 | 7, 8, 9, 10 | 22 |
| 7 | 18 | 8, 9 (+11, 12) | 3 (+14, 42) |
| 8 | 8 | 9, 10 (+11, 12) | 3 (+2, 8) |
| 9 | 0 | 10, 11 | **all 30 failed at 11** |

The window-7 and window-8 survivors were then checked at lengths 11 and 12 as
well, since two lengths is the minimum that rules out an end effect and those had
only two: **52 checks, no failures** (14 and 42 confirmations at 11 and 12 for a
window of 7; 2 and 8 for a window of 8). All 44 are in `VERIFIED_MOVES` now,
taking the listed table from 52 rules to 96 and the table with families from 64
to 116.

**And the rule that mattered was not one of the 44.** Among them,
`(0:2) -> (3:2)` via three left mutations, next to E-010's one- and two-mutation
versions, is the third member of a family whose `d`-th member needs `d`
mutations — so discovery at any bounded depth sees only an initial segment of it.
Generating the family instead gives every member: F-020, and H-008 confirmed.
That is the return on this run, more than the 44.

Cost: roughly 30 s per (pattern, embedding) at `maxSteps=3`; the re-verification
is the expensive half, since a window of 8 wants length 12 and its 58786 LNAs.
Tests H-007, and H-007 bit back.

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
