# The out-degree 3 and 4 rows with J != 0 on the n = 8 c0 walk are real Cartan failures (5 of 5), and the Cartan defect sits exactly at the J_i support; the socle reading survives

author: skeptic · round: 034 · kind: result (referee follow-up; attempted refutation failed)
thread: T5 · bears on: E-113, E-124, E-126, H-015

## Claim

On the E-126 walk (n = 8, class 0, guarded BFS, 8 456 expansions), every gate-admitted (algebra, v) row with out-degree >= 3 and J != 0 fails `mutateAtVertex(..., checkCartan=True)`: 5 of 5 (4 of out-degree 3, 1 of out-degree 4), all with parallel out-arrows (v -> 8 two or three times). They are real failures of the mutation (rewrite not congruent to R C R^T, so not tilting), not artefacts of the gate; the gate is simply blind to J != 0 here, as at out-degree 2 (E-116). Control: 199 out-degree >= 3 rows with J = 0 (41 parallel, 158 not) are all congruent. The E-124 reading J_i = Hom(S_v, e_iA) was not refuted: in all 5 failures the support of the Cartan discrepancy R C R^T - Cartan(child) is exactly {(v, i) : J_i != 0}. Not claimed: that these 5 are distinct orbits (rows 6820/7424 and 7822/7831 look like mirror pairs, unchecked), that out-degree >= 3 failures occur outside n = 8 c0, or that the rows are reachable from an LNA by a different route (they are on a key-preserving walk, so they are).

## Evidence

Walk as `experimentalist_dimji.py`, so the same rows as the round-033 referee note (expansions 6 820, 7 424, 7 822, 7 831) plus one more:

| expansion | v | out-arrows | J_i != 0 (dim) | Cartan | discrepancy support |
|---|---|---|---|---|---|
| 6 820 | 4 | 6, 8, 8 | i = 3, 5 (1, 1) | FAILS | (4,3), (4,5) |
| 7 424 | 4 | 6, 8, 8 | i = 1, 5 | FAILS | (4,1), (4,5) |
| 7 822 | 7 | 4, 8, 8 | i = 3, 5 | FAILS | (7,3), (7,5) |
| 7 831 | 7 | 4, 8, 8 | i = 1, 5 | FAILS | (7,1), (7,5) |
| 7 836 | 7 | 6, 8, 8, 8 | i = 1, 3, 5 | FAILS | (7,1), (7,3), (7,5) |

Tally by (out-degree, arrows, J, Cartan): (3, parallel, J != 0, FAILS) 4; (4, parallel, J != 0, FAILS) 1; (3, parallel, J = 0, congruent) 41; (3, simple, J = 0, congruent) 158; (4, simple, J = 0, congruent) 1. No out-degree >= 3 row with J != 0 and a congruent Cartan; no J = 0 row failing. So within this walk "J != 0 iff Cartan fails" holds at out-degree 3, 4 (205 rows; positives are the 5).

Socle attempt. J is computed as ker of x -> (x b)_b on e_iAe_v, which is the socle statement by E-124's own derivation, so agreeing with J is not independent. The Cartan discrepancy is computed from the Cartan matrices of parent and child and does not use J. Its nonzero entries are all in row v, at the columns i with J_i != 0, in each of the 5 rows (E-124 saw this at one row, E-068 step 7). One more observation, 5 rows only: the number of i with J_i != 0 equals the number of parallel arrows into 8 (2, 2, 2, 2, 3), each dim 1. Probably the count of independent "differences of parallel arrows", not tested.

Caveat on these rows: the relations are messy (repeated paths inside a relation, a relation `[2,4,8] - [2,4,8]` with different parallel arrows, three-term relations), so they are a different regime from the scalar-1 two-term W/D' rows of E-113..E-122; the E-112 circuit lemma (monomials and p - q) does not apply as stated. I did not classify Gamma_i or test W on them.

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/034/skeptic_outdeg3.py 8 0 560 9000   # 560 s; output skeptic_outdeg3_n8c0.txt
```
The cap is by time, so a loaded machine expands fewer algebras and may miss the later rows (the first appears at expansion 6 820; `max_exp` 9000 is the cap, 560 s the budget).

## Prior record

E-113: "out-degree >= 3: 3 494 rows, 0 rejects" in a 540 s prefix walk; consistent (first row is at expansion 6 820). E-126 referee note: at least 4 such rows, "run checkCartan": done here, answer is real failures. E-116/E-122 explain the out-degree 2 analogue (gate tests single paths). E-124: socle reading. Nothing in RETRACTIONS touched. New: the out-degree >= 3 failures exist, so "rejects are out-degree 1 or 2 only" (implicit in E-105/E-112/E-113 tables) is false on walks at n = 8 c0; the gate is not complete at out-degree 3, 4 either.

## Code changed

None in `quivermutation/`. New script `workshop/rounds/034/skeptic_outdeg3.py` (reuses the round-033 helpers by exec); not tested beyond the run above.

## Next

- toolsmith: does the capped walk at c1 or n = 9 c0 prefix have out-degree >= 3 rejects; the walk found its first at expansion 6 820, so use a checkpointed run.
- theorist: why the number of J_i != 0 equals the parallel multiplicity in 5 rows; Gamma_i for parallel arrows with 3-term relations (E-112 lemma needs extension); the 5 rows' orbits (mirror pairs?).
- skeptic: a socle test independent of J would need coker / Hom(N, N[-1]) (E-124 open); try a hand-built row with dim J_i >= 2.
