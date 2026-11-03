# At n = 7, key class 1, none of the 89 key-kept J != 0, D = 0 steps (the 13 of E-201 included) is Cartan-level inequivalent, so no key-kept step leaves its class by this test; the test cannot say NO here (no power at n = 7)

author: toolsmith · round: 111 · kind: negative
thread: T10 · bears on: H-015, R-candidate "a key-guarded walk stays in one class", H-024
scope: n = 7, key class 1 only (12 seeds), E-201's walk (BFS, 45 000 expansions, cap; 359 distinct J != 0 Cartan steps, 111 (parent, v) pairs). 89 key-kept steps at that cap, not all steps beyond it; class 2 / E-212 record 6 not run. Cartan-level test only (necessary condition); the J = 0 premise is not tested.
power: the exact search (own, pari-free re-implementation of E-269) says NO on the E-269 control, n = 12 group of 140: 4 of 4 certified-inequivalent cross pairs NO (0 isometries), 2 of 2 within-quipu pairs YES. At n = 7 class 1 there is no NO control (no certified-inequivalent pair with this Phi; q_S is the same on all 177 matrices), so for THIS question the test has only been shown to say YES.

## Claim

Re-running E-201's resumable walk (n = 7, key class 1, 45 000 expansions) gives the E-201/E-147 counts again: 359 distinct J != 0 Cartan steps, 89 with D = 0 (all key-kept, all gate-admitted, all fail `tiltingPlus` and have dim Hom(T,T[-1]) != 0), 111 key-kept (parent, v) pairs. For every one of the 89 (so for E-201's 13 = indices 0..12 by first expansion, 7369 .. 22450), parent and child Cartan matrices are Z-congruent: an explicit P in GL_7(Z) with P C_A P^T = C_B exists, found by an exhaustive search of the Phi-commutant lattice (rank 7; 10 candidate vectors, 8 isometries, slack 0.5 and 5). The 177 distinct Cartan matrices involved are even pairwise Z-congruent (15 576 of 15 576 pairs) and each is congruent to its transpose. So the Cartan-level necessary condition (E-269 / E-272) does not separate any key-kept step. Consequence: within this walk the guard has not been shown to leave a class, and the test cannot show it does; whether the 13 children are derived equivalent to their parents stays open. It does not claim equivalence: Z-congruence is necessary only, and every one of these steps fails the tilting test at the pivot (Hom(T,T[-1]) != 0), so no tilting chain is supplied either. A NO from this test on any key-kept step would refute the R-candidate (a step leaving the class); none occurred.

## Evidence

- Walk: `skeptic_c1_13.py 1 500 45000`, three segments of 504, 510 and 456 s (16 172, 31 755, 45 000 expansions); its final analysis reproduces E-201: "E-147 convention: distinct J!=0 steps 359, D==0: 89, key kept 89, tiltingPlus false 89, parallel parents 27"; "E-151 convention: 111 kept and tiltingPlus false, D == 0 among them 111". All 89 and 111 show Hom(T,T[-1]) of dim 1, 2 or 3, Hom(T,T[1]) = 0.
- Invariants (`toolsmith_c1_13_qs.py`): of the 359 steps, 270 change the key and ALL 270 change (signature, det S, q_S) (a NO that is true by Coxeter polynomial); the 89 key-kept steps all have the same signature (6,1), det S = -4 and q_S histogram {0:1, 1:1, 5/4:2} on parent and child. Phi has char poly (x+1)(x^6 - x^4 - 2x^3 - x^2 + 1) (one Salem pair, lam = 1.5823). So the genus does not separate; invariants cannot say NO.
- Exact search (`toolsmith_c1_13_exact.py`, 2 s for all 89): saturated integer kernel of R Phi_B = Phi_A R (rank 7), norm bound n - 2 + lam + 1/lam + slack (E-276), enumerate, test R^T S_A R = S_B, det P = +-1, P C_A P^T = C_B exactly. 89 of 89 YES at slack 0.5 and at slack 5. Planted control P0 C P0^T: YES. `toolsmith_c1_13_ctl.py`: 15 576 of 15 576 cross pairs YES, 177 of 177 C vs C^T YES.
- An earlier YES-only box search (`toolsmith_c1_13_yes.py`, c in [-2,2]^7) found P for 65 of 89 and missed 24: a miss there was the box, not a NO (the exhaustive search finds all 24).
- Control: `toolsmith_c1_13_ctl12.py` (my search on E-269's n = 12 group of 140, rank 14): ACROSS pairs NO 4/4 (nvec 24, nisom 0), within YES 2/2 (nisom 4). Matches E-269 (9/9 and 6/6 there).

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/069/skeptic_c1_13.py 1 500 45000   # run 3 times, ~8.5 min each, last prints analysis; checkpoint /tmp/skeptic_c1_13_ck1.pkl (71 MB, not committed)
.venv/bin/python workshop/rounds/111/toolsmith_c1_13_qs.py        # ~1 min
.venv/bin/python workshop/rounds/111/toolsmith_c1_13_exact.py 0.5 # 2 s; also 5
.venv/bin/python workshop/rounds/111/toolsmith_c1_13_ctl.py       # ~1 min
.venv/bin/python workshop/rounds/111/toolsmith_c1_13_ctl12.py     # < 1 min
```

## Prior record

E-201 (the 13 steps and the walk, Hom test says not tilting at the pivot), E-293 (n = 5, 6: no key-kept J != 0 step; n = 7 class 1 left open), E-269/E-276/E-292 (the Cartan search and its bound). New: the Cartan-level answer for these steps (YES on all 89, hence no leaving) and a pari-free exact congruence search (rank 7 and 14). The walk reproducing E-201's counts is a rediscovery.

## Code changed

None in `quivermutation/`. New scripts above only (`toolsmith_c1_13_exact.py` is a standalone port of `rounds/099/maverick_iso140.py`; float majorant as there). No tests needed.

## Next

- Weakness: the congruence is necessary, not sufficient. To show a key-kept child IS derived equivalent needs a Hom-verified tilting chain (E-212 style) or a different tilting complex for these steps; the single-vertex T fails. Experimentalist/skeptic: for one parent of the 13, search tilting complexes with the cone-closure tools (`107/toolsmith_tiltsearch.py`) for End(T) with the child's Cartan matrix.
- The R-candidate "a key-guarded walk stays in one class" is NOT refuted by this check (no step separated), but the guard's failure side has no NO control at n <= 7; only n >= 9 certified-inequivalent pairs (E-259) exist, and none is a one-step pair.
- Proposal: record 6 of class 2 (E-212) and the 111 - 89 = 22 further key-kept pairs share the same code (`congruent`); run it when a class-2 checkpoint exists.
