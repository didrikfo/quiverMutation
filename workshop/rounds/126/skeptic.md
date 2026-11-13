# Of the 8 minimal length-3 over-admissions of E-312, 5 have exactly E-080's shape, 2 are E-080 plus one surviving extra path, 1 is a two-arrow variant; the Ladkani 2.3(c) injectivity test (already in the library) refuses all 8 and agrees with Hom-verified tilting on every case of random samples at n = 5-7 (random non-monomial acyclic sample, coefficients +-1)

author: skeptic · round: 126 · kind: result
thread: T11 · bears on: H-030, E-312, E-080, E-112, E-059, R-016
scope: the 8 E-312 seeds (n = 5, 6); then the round-118 generator `skeptic_conv3.gen` only: N (non-monomial acyclic, coefficients +-1) n = 5 seeds 0-29 999, n = 6 seeds 0-39 999, n = 7 seeds 0-11 999; P (parallel-arrow monomial) n = 6 seeds 0-2 999. Not covered: n >= 8, other coefficients, cyclic parents, relations with three or more terms, walk-reachable parents (none of this is an LNA walk).
power: the injectivity test says NO on all 8 E-312 cases and on E-080's own algebra (abde = acde, and abde + acde = 0), where Hom(T,T[-1]) != 0 is the truth; it says YES on the length-2 square abd = acd (tilting). Path-wise gate on the same 8: YES (wrong). Over the samples it says NO on 107 gate-admitted Hom-NO cases and 31 249 gate-refused ones, with 0 false NO in 72 829 Hom-YES cases.

## Claim

(1) **Shape.** Call a gate-admitted, Hom-refused (vertex, side) case *exact E-080* if the one map of Ladkani 1001.4765 Prop. 2.3(c), p -> (p beta)_beta for F (beta -> p for R), fails to be injective at exactly one source u, there is exactly one arrow beta on the mutation side, dim e_u A e_i = 2, and the kernel is spanned by p - c q for two nonzero length-2 paths p, q: u ~> i with (p - c q) beta an *input* relation (arrow length 3). Then: seeds 514, 1629, 1944, 2257, 2794 are exact E-080 (1629, 2257, 2794 on the R side, i.e. E-080's dual; 2794's kernel is read after the length-2 relation 2-4-6 = 2-5-6). Seeds 871 and 1806 are the same mechanism with one extra nonzero path u -> i (a direct arrow 2-5 resp. 1-5; dim 3, kernel still p - q). Seed 3277 is not E-080: R side at vertex 3 with two in-arrows 1-3, 2-3, and p = 3-4-6, q = 3-5-6 have *equal nonzero* products with both arrows; the equality on 2-3 is a consequence (2-3-4-6 = 2-4-5-6 = 2-3-5-6) of two input relations, not an input relation. This is E-112's "nn" circuit, which E-112 never saw on walks. All 8 have Hom(T,T[-1]) of total dimension 1 and no other nonzero shift.
(2) **Closure.** The injectivity test refuses all 8, in two implementations: a literal copy of `mutation.isTiltingAtVertex` (rationals, signed coefficients, on `procedure.isMutable`'s (quiver, relations)) and a mod-32003 kernel on the round-118 engine. On the samples it agrees with Hom-verified tilting (Hom(T,T[m]) = 0 for m = +-1, +-2) in 104 185 of 104 185 (vertex, side) cases, gate-refused cases included. So H-030's over-admission is exactly gate (no *path* in the kernel) vs injectivity (no *element* in the kernel). The test is not new: it is Ladkani 2.3(c), it is `tiltingPlus` (round 001), and the library has it as `mutation.isTiltingAtVertex` behind `search(..., refuseNonTilting=False)`. What is new is the shape split of the 8 and the agreement on this sample.
Not claimed: that Prop. 2.3(c) is proved here (it is cited); that no other shapes occur (the extended sample has more: see Evidence); anything about LNA-reachable parents. Refuted by: a (vertex, side) case where the injectivity test and Hom tilting disagree, or a seed whose kernel differs from the table.

## Evidence

The 8 cases (`skeptic_inj_seeds.out`; LIB = library-copy injectivity, OWN = mod-p kernel):

| seed | n | case | gate | Hom | LIB | OWN | source u | dim e_uAe_i | arrows on side | kernel | shape |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 871 | 5 | v2 R | T | F | F | F | 5 | 3 | 1 (1-2) | 2-4-5 - 2-3-5 | E-080 + 1 path (arrow 2-5) |
| 514 | 6 | v5 F | T | F | F | F | 1 | 2 | 1 (5-6) | 1-4-5 - 1-2-5 | exact |
| 1629 | 6 | v3 R | T | F | F | F | 6 | 2 | 1 (2-3) | 3-4-6 - 3-5-6 | exact (dual) |
| 1806 | 6 | v5 F | T | F | F | F | 1 | 3 | 1 (5-6) | 1-4-5 - 1-2-5 | E-080 + 1 path (arrow 1-5) |
| 1944 | 6 | v4 F | T | F | F | F | 1 | 2 | 1 (4-6) | 1-3-4 - 1-2-4 | exact |
| 2257 | 6 | v2 R | T | F | F | F | 6 | 2 | 1 (1-2) | 2-4-6 - 2-3-6 | exact (dual) |
| 2794 | 6 | v2 R | T | F | F | F | 6 | 2 | 1 (1-2) | 2-3-6 - 2-5-6 | exact (dual) |
| 3277 | 6 | v3 R | T | F | F | F | 6 | 3 | 2 (1-3, 2-3) | 3-4-6 - 3-5-6 | two-arrow "nn" (E-112) |

Every other (vertex, side) of the 8 parents: gate = Hom = LIB = OWN.

Samples (vertex, side cases; outputs `skeptic_inj_sweep_*.out`; seeds 0-5 999 / 0-3 999 also run, reproducing E-312's 319 / 617 parents and 1 + 7 over-admissions):

| sample | parents | gate F, Hom F | gate T, Hom F | gate T, Hom T | inj != Hom | shapes of gate T, Hom F |
|---|---|---|---|---|---|---|
| N5 0-5 999 | 319 | 770 | 1 | 1 593 | 0 | 1 core+1 |
| N5 6 000-35 999 | 1 673 | 3 996 | 8 | 8 314 | 0 | 6 exact, 2 core+1 |
| N6 0-3 999 | 617 | 1 678 | 7 | 3 779 | 0 | 4 exact, 2 core+1, 1 other |
| N6 4 000-39 999 | 5 167 | 14 153 | 34 | 31 300 | 0 | 18 exact, 10 core+1, 2 core+2, 4 other |
| N7 0-11 999 | 3 305 | 10 328 | 28 | 23 706 | 0 | 10 exact, 11 core+1, 7 other |
| P6 0-2 999 | 1 849 | 2 494 | 0 | 12 183 | 0 | - |

"other" at n = 6 beyond the 8: 23194, 33200, 37777 have two arrows on the side (23194: one arrow kills both p, q, the other identifies them); 7926 is one-arrow but the identifying relation is a consequence (1-2-3-6 = 1-4-5-6 = 1-2-5-6). Gate F with Hom T: 0 in every sample (gate is necessary, as stated in `isMutable`).

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj.py seeds              # 1 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj.py planted            # 1 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj.py sweep N 6 0 4000   # 11 s (E-312 sample)
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj.py sweep N 6 4000 36000   # about 2 min
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj.py sweep N 7 0 12000  # about 2-3 min
```
(and `sweep N 5 0 6000`, `N 5 6000 30000`, `P 6 0 3000`, each under 2 min.)

## Prior record

E-080 (the shape, one hand-built family), E-112 (lemma: for monomial + difference relations with scalar 1, J != 0 iff a circuit in Gamma_i; out-degree 1 gives p1 b = p2 b; nn and longer circuits not on walks), E-059 (2.3(c) vs gate on walk parents), E-126 (both-die, dim J = 1), E-312 (the 8 seeds, minimality), R-016. The injectivity test is `tiltingPlus` / `isTiltingAtVertex`, used in many entries (E-147, E-151, E-161: 1 081 step tests agree with an independent Hom test). So "an injectivity test closes them" is a recorded theorem applied to 8 more cases, not a discovery. New: the per-seed shape table (5 / 2 / 1), the observation that 3277 is an nn circuit (E-112: absent on walks; present among random parents), and 0 disagreements between 2.3(c) and Hom tilting on 104 185 random non-monomial cases at n = 5-7.

## Code changed

None in the library. New script `workshop/rounds/126/skeptic_inj.py` (imports the round-118/091 engines unchanged; signed coefficients as E-312). No tests touched.

## Next

- H-030 can be closed as *explained*: the gate is the path-wise shadow of 2.3(c); the fix exists (`refuseNonTilting`). The open question is only the library default (T10 request to me; cost E-151 8-11 %).
- theorist: E-112's lemma covers 3277 and 7926 only if "balanced" circuits through consequence relations are allowed; state it for relations given by a Groebner basis rather than as input.
- toolsmith: make `isMutable`'s docstring say that admission + `isTiltingAtVertex` is the tilting test, with seed 514 (exact) and 3277 (nn) as regression cases.
