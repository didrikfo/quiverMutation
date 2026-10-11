# At n = 7 (classes 1, 2), End(T) is isomorphic to the next algebra on all 25 failing n = 7 steps (8 parallel-arrow ones decided) and on 405 path edges; the comparison did not discriminate at any of the 25 failing steps (Hom(T,T[-1]) != 0 at all 25)

author: toolsmith · round: 057 · kind: result (with a negative half)
thread: T10 · bears on: E-157, E-160, E-161, E-167, E-168, E-169, H-015 (indirectly)
scope: n = 7 classes 1, 2 (E-151 walk rebuilt at 20 000 expansions; keys of all 46 logged children/parents reproduced) and n = 10 (the 5 E-169 paths); label-preserving iso over the algebraic closure; Hom_K exact over Q; generation of K^b(proj) by T not re-tested here (E-168 for loop-free steps); not covered: n = 8, the 25 E-157 / E-160 paths for J != 0 (none exist), class-2 negative control.

## Response to referee
Verdict received: minor revision. New script `toolsmith_endt2_hom.py` (outputs `toolsmith_endt2_hom_c1.txt`, `_c2.txt`; a few seconds each).
1. **Done.** The skeptic's `stepTest` (rounds/050/skeptic_tilt.py) gives Hom(T,T[-1]) = 1 (dimension sum) and Hom(T,T[1]) = 0 at **16 of 16 c1 and 9 of 9 c2** failing steps, so "25 of 25 J != 0" is now computed, not asserted (c2 was not in E-168's check).
2. **Done.** Title and Claim now say End(T) "did not discriminate at 25 of 25 failing steps"; the lemma "End(T) = mutation algebra at every step" is a conjecture on n = 7 and n = 10 (one start), not a theorem.
3. **Done.** Claim now states that the dimensions of K Q_c/I_c (`child_info`) and the completeness of `crels` as a generating set of I_c are inputs of the test, not outputs; a dropped relation in `reduction` with equal dims would pass.
4. Next round (proposal kept): wrong-algebra control with equal dims, arrows and a parallel pair (one relation altered at the same Cartan matrix).

## Claim
(1) The 8 parallel-arrow failing c1 steps left undecided by E-167 are decided: End(T) (T = P_v -> sum P_h, S = -1) is isomorphic, label-preserving, to the mutated algebra at each of them (steps 0, 1, 2, 3, 8, 10, 11, 12; multiplicity 2, two parallel pairs in step 12; path edges reach multiplicity 3 once). With the 8 of E-167, **16 of 16 c1 and 9 of 9 c2 failing J != 0 steps are decided, all iso**. (2) On the E-157 (40 paths), E-160 (3) and E-163 (path13) paths of n = 7: **370 of 370 edges iso** (c1 214, c2 156; 56 c1 edges have parallel arrows, 0 undecided); the five E-169 n = 10 paths: **35 of 35 iso** (no parallel arrows). Undecided: 0; non-iso: 0. The refuting event would be one edge with the Groebner ideal containing 1 or differing arrow counts; none occurred.
Inputs, not outputs, of the test: the dimensions of K Q_c/I_c (`child_info`) and that `crels` generates the whole ideal I_c; an incomplete `crels` with equal dims would pass undetected. It does **not** say the J = 0 premise holds: because End(T) is also iso to the child at every J != 0 step, this comparison did not discriminate at any of the 25 failing steps, all of which have Hom(T,T[-1]) != 0 (computed, 16 + 9; E-167 negative half now 25 of 25, not 8 of 16). Class 2 and n = 10 added only agreement of two steps of the same kind.

## Evidence
Method (`toolsmith_endt2.py`, `symcheck2`): for each pair (i,j) with m arrows of c the arrows are sent to M * (lifts of a basis of rad/rad^2) + rad^2 corrections, M an unknown m x m matrix; every relation of c must vanish in End(T); one Rabinowitsch variable per block enforces det M != 0 (z * product of all dets was too slow: one 3 x 3 + 2 x 2 block hung sympy for 10 min); Groebner basis != [1] => a solution exists over the closure; with equal dims and arrows generating (E-167 check), K Q_c/I_c -> End(T) is an isomorphism. For m = 1 everywhere this is the round-053 test (same outputs on all 13 path13 edges and the 16 - 8 old decided steps).

| set | steps/edges | iso | parallel-arrow | undecided |
|---|---|---|---|---|
| failing J != 0, c1 | 16 | 16 | 8 | 0 |
| failing J != 0, c2 | 9 | 9 | 0 | 0 |
| E-157 paths c1 (23) + E-160 (2) + E161 path13 (1) | 214 edges | 214 | 56 | 0 |
| E-157 + E-160 paths c2 (17 + 1) | 156 edges | 156 | 0 | 0 |
| E-169, 5 paths, n = 10 | 35 edges | 35 | 0 | 0 |

Power control of the matrix-valued test (`toolsmith_endt2_ctrl.py`, the 8 parallel steps, 9 multi-term relations): replacing a binomial by its first term is rejected 9 of 9; doubling a coefficient is absorbed 9 of 9 (expected: the scalar blocks absorb it, as in E-167). Weakness: no wrong-algebra control with equal dims, arrows and a parallel pair; the relation-level power is the monomialisation control above. The E-167 caveat "wrong-algebra control only at the dims filter" stands.
Max Hom dimension on a path edge is 4 (c1); n = 10 edges all have dim <= 1 (so that test is the weakest: 35 edges with maxdim 1 are decided by dims and arrows plus a few scalar equations).
Determinism check: the rebuilt pickles reproduce the 28 (c1) and 18 (c2) child/parent keys printed in `rounds/049/toolsmith_paths_logs.txt` (`toolsmith_endt2_keys.py`, 0 mismatches).
Outputs: `toolsmith_endt2_{fail_c1,fail_c2,paths_c1,paths_c2,n10,ctrl}_out.txt` (total < 100 KB).

## Reproduction
```
for c in 1 2; do DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 $c 20000 /tmp/tsm/c$c.pkl 100; done   # c2 once (480 s), c1 twice (about 520 s + 120 s)
CLS=1 .venv/bin/python workshop/rounds/057/toolsmith_endt2_keys.py                       # also CLS=2; keys match
CLS=1 .venv/bin/python workshop/rounds/057/toolsmith_endt2_run.py plan                   # 26 c1 paths, 214 edges (c2: 18 paths, 156)
CLS=1 timeout 10m .venv/bin/python workshop/rounds/057/toolsmith_endt2_run.py fail       # 16 steps, < 5 s; CLS=2: 9 steps
CLS=1 timeout 10m .venv/bin/python workshop/rounds/057/toolsmith_endt2_run.py paths 0:40 # whole c1 list, about 1 minute; CLS=2 likewise
timeout 10m .venv/bin/python workshop/rounds/057/toolsmith_endt2_n10.py run 0:5          # 35 edges, < 1 minute
CLS=1 timeout 10m .venv/bin/python workshop/rounds/057/toolsmith_endt2_ctrl.py           # perturbation control
```
`--plan` first: `plan` lists paths and edge counts without computing End(T). No job exceeded 2 minutes, so no sharding was needed.

## Prior record
E-167 (13 edges, 8 of 16 failing steps decided, 8 undecided) is extended, not contradicted: its decided cases and path13 give the same verdicts here. E-161/E-166 are the Cartan-level Hom test on the same edges; E-168 closes generation per loop-free step. E-169 recorded the n = 10 paths without an End(T) test; this adds one. Nothing in RETRACTIONS.md bears on it (grepped "End(T)", "parallel"). The earlier note "sympy groebner with z * prod" is a performance trap, not a result.

## Code changed
None in `quivermutation/`. New scripts only (`toolsmith_endt2*.py`, which reuse the round-053 `TiltEnd`). No tests run (no library file touched).

## Next
- Discriminating the premise needs the extension group, not End(T): the useful check is the Hom(T,T[-1]) = 0 test at the n = 10 E-169 edges (E-169 said the gate never fails there; the skeptic's `stepTest` is n-independent). I did not run it.
- theorist: the observation that End(T) = mutation algebra at J != 0 steps too (25 of 25) is a lemma candidate (mutation algebra = End of the 2-term complex always; J measures only Hom(T,T[-1])). If proved, the E-167 comparison can be dropped and T10 reduces to Hom(T,T[-1]) = 0 per step.
- skeptic: replay `symcheck2` on two parallel-arrow edges by hand (e.g. fail c1 step 0, matrix 2 x 2 at pair (4,5)).
- toolsmith (later): `DEFAULT_CAP` 5040 and docstring rewords remain open.
