# The Cartan congruence of E-085 and `tiltingPlus` are the same condition on every gate-admitted step tested: congruence fails exactly where the one map of AI 2.32(b) / Ladkani 2.3(c) has a kernel, and the discrepancy equals minus its dimension

author: scholar · round: 018 · kind: result (with a short derivation, one step conditional)
thread: T5 · bears on: H-015, E-066, E-078, E-085, E-089

## Claim

Write g_i : e_iAe_k -> (+)_beta e_iAe_{head beta}, p |-> (p beta)_beta, the map of AI 2.32(b) / Ladkani 2.3(c)
(`tiltingPlus` = "g_i injective for all i != k"). Let X = R C R^T (Ladkani Prop 3.6 / Lemma 3.5, the Euler form of the two-term
complex T^+_k) and Y = Cartan of the child the rewrite produces. Then, on all steps tested, (i) tiltingPlus True <=> X = Y;
(ii) when X != Y the difference X - Y is supported in row k, off the diagonal, with (X - Y)_{k,i} = -dim ker g_i exactly.
Mechanism: (RCR^T)_{k,i} = chi(P_i, T^+_k) = dim coker g_i - dim ker g_i, while the rewrite's child has dim coker g_i
(its e_iA'e_k is the cokernel of g_i). So the Cartan test is the **row-k off-diagonal part of "g_i injective", read through the
rewrite**; it is not an independent second criterion. Not claimed: that congruence implies tilting for an arbitrary
algebra (X = chi, which counts Hom^0 - Hom^{-1} - Hom^1, can in principle cancel); that the child is derived equivalent
(that rests on Ladkani's theorem); that the one-map identity AI = Ladkani = `tiltingPlus` is derived (still the authors' word).

## Evidence

Library as after E-089 (full-reduction `reduceAgainstPivots`). Crosstab (tilt, congruence) at every gate-admitted step
of a guarded BFS from the LNAs and relation duals, steps beyond the guard also tested (every child, guard or not):

| set | algebras | steps tilt+cong | steps NOT tilt, NOT cong | xor | NOT-cong with diff = -dim ker (row k only) |
|---|---|---|---|---|---|
| E-078 family n = 5..7 (18 algebras, all vertices) | 18 | 75 | 6 | 0 | 6 of 6 |
| n = 5, both key classes, closed | 11 700 | 30 300 | 0 | 0 | - |
| n = 6 class 0, stopped at 420 s (not closed) | 32 450 | 65 914 | 696 | 0 | 696 of 696 |
| n = 7 class 0, stopped at 420 s (not closed) | 22 016 | 39 901 | 111 | 0 | 111 of 111 |

All 807 non-tilting steps are guard-refused (consistent with E-084). The 10 n = 8 class 2 "key moved" steps of E-085 (a rewrite
defect, congruence failed while `tiltingPlus` was True) are the only recorded disagreement; with the fixed rewrite none appears
in these runs, and the n = 8 class 2 walk itself is not re-run here (E-090).
E-078 at n = 5, vertex d: X - Y has a single nonzero entry, -1 at (d, a), and ker g_a = <abd - acd> has dimension 1. The
`d`-row of X is (-1,0,0,1,1), of Y (0,0,0,1,1): the rewrite cannot represent a negative Hom, so it drops to 0.

Derivation sketch. T^+_k = (P_k -> (+)P_beta). Apply Hom(P_i, -): H^{-1} = ker g_i, H^0 = coker g_i, so
chi(P_i, T^+_k) = dim coker g_i - dim ker g_i. For a tilting complex all of Hom(T,T[m]), m != 0, vanish and the Cartan
matrix of End T is chi; so tilting => X = Y (Lemma 3.5). Conversely, with the rewrite giving dim coker g_i in position (k,i),
X = Y at (k,i) iff ker g_i = 0. **Least sure step:** that the rewrite's (k,i) entry is dim coker g_i in every case; it is read
off the data above (807 of 807), not from the code of `procedure.mutateAtVertex` step 7.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/018/scholar_cartan_vs_tilt.py --e078          # 3 s
timeout 10m .venv/bin/python workshop/rounds/018/scholar_cartan_vs_tilt.py 5 --all         # 86 s
timeout 10m .venv/bin/python workshop/rounds/018/scholar_cartan_vs_tilt.py 6 --class 0 --budget-sec 420   # 7 min
timeout 10m .venv/bin/python workshop/rounds/018/scholar_cartan_vs_tilt.py 7 --class 0 --budget-sec 420   # 7 min
.venv/bin/python workshop/rounds/018/scholar_e078_diff.py                                  # 2 s
```
Outputs: `workshop/rounds/018/scholar_cartan_n6_c0.txt`, `scholar_cartan_n7_c0.txt`.

## Prior record

E-055 ran both tests on 61 718 steps (n = 6, 7) with 0 failures but as two separate checks, never as one question; E-078 and
E-085 show both failing together (17 parents) and E-085 leaves "one criterion?" open (STATE T5). Lit note
`literature/1001.4765` already says Prop 3.6 is the specialisation of Lemma 3.5. New here: the identification of the
discrepancy with dim ker g_i, the row-k localisation, and the 807 non-tilting steps (E-084 had 0 Cartan data for them beyond
the 11 replays). Consequence: E-085's "Cartan congruence agrees with `tiltingPlus`" is not two pieces of evidence for
non-tilting; and on the tilting side congruence is a check on the **rewrite** (it caught the E-085 defect), which `tiltingPlus` cannot do.

## Code changed

None in the library. New scripts `scholar_cartan_vs_tilt.py`, `scholar_e078_diff.py` (no tests: analysis only).

## Next

- Toolsmith: Cartan congruence as an assertion inside `mutateAtVertex` (cheap, catches rewrite defects like E-085); keeps `tiltingPlus` as the gate.
- Theorist: confirm "rewrite's (k,i) entry = dim coker g_i" from step 7, and whether diagonal / column-k entries can ever differ.
- Anyone: an algebra with gate True where ker g_i != 0 yet X = Y would refute the equivalence; none seen (807 cases).
