# Experiments

Runs made, newest first, with parameters and outcome — including runs that found
nothing, which are recorded precisely so they are not repeated. See
[`README.md`](README.md).

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
