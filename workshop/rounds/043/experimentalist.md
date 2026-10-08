# At n = 8 no J != 0 step has the E-143 shape (H1, H2): all 32 distinct steps found have |out v| = 2 (or |supp J| = 2), so s and c_2 do not exist there, and in class 1 Q's lowest term is x^3 in 5 of 7 steps

author: experimentalist · round: 043 · kind: negative
thread: T5 · bears on: E-143, E-141, E-140, E-138

## Claim

Agenda item 1, the n = 8 tally for E-143 (guarded class walk as E-143, 520 s cap per cell, 4 cores shared by 6 jobs).
(a) H1/H2 steps: n = 8 classes 0, 1, 2 have **0**. Every J != 0 step found at n = 8 (c0: 20, c1: 7, c2: 0) is off-shape: c0 20/20 have |out v| = 2, |supp J| = 1 (so C'e_v = e_v holds but u has two arrows out); c1 has 5 with |out v| = 2 and 2 with |supp J| = 2 (J = (1,1)). So the E-143 reduction (P1, P2) is vacuous at n = 8: s, c_2 are undefined there. E-143's shape (|out v| = 1) is a property of n = 6, 7, not of the walk.
(b) The lowest-term law of E-141 ("lowest term x^2, coefficient 1") holds at n = 8 c0 (20/20) but **fails in c1**: Q has lowest term x^3 in 5 of 7 steps (coefficient 1 in 4, 2 in 1) and x^2 coeff 1 in 2. It is the first record of a lowest-degree-3 Q. It does not give a key-preserving step (Q != 0 in all 7; consistent with E-140). So "lowest term x^2" is no longer a law across n; it is n = 6, 7 (and n = 8 c0).
(c) Same script on n = 6, 7 (time-capped) reproduces E-143's shape: n = 6 c0 145 of 146 steps H1&H2 (104 distinct (Z,w,i)), P1 formula = library Q in 145/145, Q lowest x^2 coeff 1 in 146/146, F of order 8; s = 1 in 56 (c_2 = 0 in all 56), s = -2 (= 6 mod 8) in 89 with c_2 = 1. n = 7 c0 49/49 H1&H2 (38 distinct), F order 5, s in {-4, 1, 6} contains 1 (mod 5), c_2 = 0 in 49/49. n = 7 c1: 10/10 H1&H2, F order 12, s = 1 in 2 and **no s at all in 8** (the window |s| <= 8 covers every residue mod 12, so e_i is not in the F-orbit of e_w there; E-143's "never absent" held only for n = 6, 7 c0). No c_2 != 0 with s = 1 anywhere (0 of 56 + 49 + 2).

## Evidence

Cell = (n, class, seeds); steps = gate-admitted with J = perI != 0; distinct = distinct (C_A, v).

| cell | seeds | expanded / depth | J != 0 steps (distinct) | H1&H2 | off-shape shapes (|out v|, |supp J|, dim J) | Q lowest (deg, coeff) |
|---|---|---|---|---|---|---|
| n6 c0 | 2 | 7342 / 11 | 146 (146) | 145 | 1 with u = e_a + e_b - e_i | (2,1) x146 |
| n7 c0 | 2 | 4517 / 8 | 49 (49) | 49 | 0 | (2,1) x49 |
| n7 c1 | 8 | 6490 / 8 | 10 (10) | 10 | 0 | (2,1) x10 |
| n8 c0 | 2 | 2448 / 9 | 20 (20) | 0 | 20 x (2,1,(1,)) | (2,1) x20 |
| n8 c1 | 8 | 3895 / 7 | 7 (7) | 0 | 5 x (2,1,(1,)), 2 x (1,2,(1,1)) | (2,1) x2, (3,1) x4, (3,2) x1 |
| n8 c2 | 18 | 4702 / 7 | 0 | 0 | none | none |

All six runs hit the 520 s cap (counts depend on load, a few percent; n = 8 walks reached only depth 7-9 with 2400-4700 algebras expanded). A cap is not a verdict: n = 8 has 0 H1/H2 steps in the explored prefix, not in the class. Honest comparison with E-140 (n = 8 c0 and c1 had J != 0 steps, c2 none; here the same cells, same emptiness for c2) and E-138 (192 J != 0 steps at n = 8 c0 with out(i) = 3): E-138 got far more steps (192 vs 20) with a longer/guard-off walk; my 20 are a shallower sample. Nothing here contradicts E-138's "out(i) = 3": |out v| = 2 here, a different quantity.
The "n = 7 c1 s absent" row: 8 steps, Z order-12 F; reported from `s set (|s|<=8)` = () in the output file.

## Reproduction

```
.venv/bin/python -u workshop/rounds/043/experimentalist_tally.py 8 1 1 --plan     # classes: sizes 2, 8, 18, 20, 26, 52 (n = 8)
timeout 10m .venv/bin/python -u workshop/rounds/043/experimentalist_tally.py 8 1 520 > workshop/rounds/043/experimentalist_tally_n8c1.txt   # 520 s; cells c0, c1, c2 and n = 7 c0 c1, n = 6 c0 the same way
```
Outputs: `workshop/rounds/043/experimentalist_tally_n{6c0,7c0,7c1,8c0,8c1,8c2}.txt`. Six jobs ran at once on 4 cores.

## Prior record

E-143 (n = 6, 7 only; asked for this tally), E-141 (lowest term x^2, n = 6, 7), E-140 (J != 0 cells, key guard off), E-138 (n = 8 c0, out(i) = 3). Grep of EXPERIMENTS/FINDINGS/RETRACTIONS/HYPOTHESES for "lowest term"/"lowest degree": only E-141, E-143. A lowest-degree-3 Q and the |out v| = 2 shape at n = 8 are not recorded. Not retracted anywhere.

## Code changed

New: `workshop/rounds/043/experimentalist_tally.py` (dump + E-143 reduce in one pass, `--plan`, `off DEPTH` option untested). No library file touched; no tests run.

## Next

- Theorist: the E-143 reduction for |out v| = 2 (u = e_w1 + e_w2 - e_i): Schur block is the same with a rank-2 row; does Q's lowest term x^2 become coeff 1 iff a pair of moments vanish? Explain the x^3 cases in n = 8 c1.
- Skeptic: a hand check of one of the 4 (3,1) steps (print them: `ex` list, not saved by me); is the c1 child still derived-equivalent? (Q lowest term 3 is the first Q_2 = 0.)
- Experimentalist next round: n = 8 c0/c1 guard-off depth 6 to raise the 20/7 counts (needs `off` mode tested, > 10 min: OVERNIGHT.md candidate), n = 9 c0 shape (|out v|?).
