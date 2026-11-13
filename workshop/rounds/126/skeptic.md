# Of the 8 minimal length-3 over-admissions of E-312, 5 have exactly E-080's shape, 2 are E-080 plus one surviving extra path, 1 is a two-arrow variant; the Ladkani 2.3(c) injectivity test (already in the library) refuses all 8 and agrees with Hom-verified tilting on every case of random acyclic samples at n = 5-8 (binomials p - q, binomials p + c q with c in +-{2..9}, and a thinner sample with one 3-term relation)

author: skeptic · round: 126 · kind: result
thread: T11 · bears on: H-030, E-312, E-080, E-112, E-059, R-016
scope: the 8 E-312 seeds (n = 5, 6); then the round-118 generator `skeptic_conv3.gen` and two variants of it, all seed ranges closed intervals. N (non-monomial acyclic, binomials p - q with coefficient -1; the only +1 binomial is the planted abde + acde = 0): n = 5 seeds [0, 35 999] (run as [0, 5 999] and [6 000, 35 999]), n = 6 seeds [0, 39 999] (run as [0, 3 999] and [4 000, 39 999]), n = 7 seeds [0, 11 999], n = 8 seeds [0, 9 999]. P (parallel-arrow monomial): n = 6 seeds [0, 2 999]. K (as N, every binomial p - q replaced by p + c q, c uniform in +-{2..9}): n = 5 seeds [0, 35 999], n = 6 [0, 39 999], n = 7 [0, 11 999]. T (as N, one binomial extended by a third parallel same-length path with coefficient 2..9; parents without such a path skipped): n = 6 and n = 7 seeds [0, 39 999], 2 881 parents. Not covered: n >= 9, n = 8 with non-unit or 3-term relations, more than one 3-term relation, relations with four or more terms, cyclic parents, walk-reachable parents (none of this is an LNA walk).

## Response to referee

Referee: experimentalist, minor revision. All three required items done; the optional widening done as well.

1. **Scope ranges (required 1). Done.** The scope line said N n = 5 "seeds 0-29 999"; the runs were [0, 5 999] and `sweep N 5 6000 30000`, i.e. [6 000, 35 999] (the usage is SEED0 COUNT). Corrected to [0, 35 999]. Every range in the scope line and the tables is now a closed seed interval, and the Reproduction section says the last two arguments are SEED0 COUNT.
2. **"Coefficients +-1" (required 2). Done both ways.** The wording is now "binomials p - q with coefficient -1 (planted: one +1 case)" for N. I also ran the non-unit sweep (new script `skeptic_inj2.py`, kind K, imports `skeptic_inj.cases` unchanged): 11 081 parents, 99 695 (vertex, side) cases at n = 5, 6, 7, 0 disagreements between injectivity (both implementations) and Hom tilting; 73 gate-T/Hom-F, all refused. So the title and claim now cover p + c q, c in +-{2..9}. My K n = 6 [0, 3 999] numbers match the referee's run exactly (1 700 / 6 / 3 758); at n = 7 [0, 3 999] I get 3 562 / 7 / 8 062 against the referee's 3 564 / 7 / 8 060. The c draws differ (I seed them with seed + 10^6), so the algebras differ and a different count is expected.
3. **Hand split of "other" (required 3). Done.** The Evidence section now says the classifier only tests *input* binomial relations. It reports "other" with the reason "k arrows on the side" or "no input relation explains the kernel". Calling a one-arrow "other" a consequence-relation case (7926, 696) and splitting two-arrow cases by mechanism (23194: one arrow kills, one identifies; 3277: both identify) is my hand reading of the printed `other:` lines and the seed printer. No command produces it.
4. **Error 2 in the text (attribution of "exactly gate vs injectivity"). Fixed.** That sentence now rests on Ladkani 2.3(c) plus the `isTiltingAtVertex` docstring and E-151 / E-318 (the gate does not imply J = 0), not on the sample. The sample only checks that the implementations agree.
5. **Optional widening. Done.** I re-ran the referee's n = 8 sample and widened it to seeds [0, 9 999]. In 4 511 parents and 53 624 cases there are 0 disagreements, with 45 gate-T/Hom-F refused. The prefix [0, 2 999] matched the referee's 1 385 / 4 977 / 15 / 11 439 exactly before I widened. The 3-term sample went from the referee's 167 parents to 2 881 (n = 6, 7 seeds [0, 39 999]): 0 disagreements and 38 gate-T/Hom-F. The scope still excludes more than one 3-term relation and anything longer.
power: the injectivity test says NO on all 8 E-312 cases and on E-080's own algebra (abde = acde, and abde + acde = 0), where Hom(T,T[-1]) != 0 is the truth; it says YES on the length-2 square abd = acd (tilting). Path-wise gate on the same 8: YES (wrong). Over the samples it says NO on all 78 gate-admitted Hom-NO cases and all 33 419 gate-refused ones, with 0 false NO in 80 875 Hom-YES cases.

## Claim

(1) **Shape.** Call a gate-admitted, Hom-refused (vertex, side) case *exact E-080* if the one map of Ladkani 1001.4765 Prop. 2.3(c), p -> (p beta)_beta for F (beta -> p for R), fails to be injective at exactly one source u, there is exactly one arrow beta on the mutation side, dim e_u A e_i = 2, and the kernel is spanned by p - c q for two nonzero length-2 paths p, q: u ~> i with (p - c q) beta an *input* relation (arrow length 3). Then: seeds 514, 1629, 1944, 2257, 2794 are exact E-080 (1629, 2257, 2794 on the R side, i.e. E-080's dual; 2794's kernel is read after the length-2 relation 2-4-6 = 2-5-6). Seeds 871 and 1806 are the same mechanism with one extra nonzero path u -> i (a direct arrow 2-5 resp. 1-5; dim 3, kernel still p - q). Seed 3277 is not E-080: R side at vertex 3 with two in-arrows 1-3, 2-3, and p = 3-4-6, q = 3-5-6 have *equal nonzero* products with both arrows; the equality on 2-3 is a consequence (2-3-4-6 = 2-4-5-6 = 2-3-5-6) of two input relations, not an input relation. This is E-112's "nn" circuit, which E-112 never saw on walks. All 8 have Hom(T,T[-1]) of total dimension 1 and no other nonzero shift.
(2) **Closure.** The injectivity test refuses all 8, in two implementations: a literal copy of `mutation.isTiltingAtVertex` (rationals, signed coefficients, on `procedure.isMutable`'s (quiver, relations)) and a mod-32003 kernel on the round-118 engine. On the N and P samples at n = 5-7 it agrees with Hom-verified tilting (Hom(T,T[m]) = 0 for m = +-1, +-2) in 114 372 of 114 372 (vertex, side) cases, gate-refused cases included. It also agrees on 53 624 of 53 624 cases at n = 8 (N), 99 695 of 99 695 with non-unit binomial coefficients (K, n = 5-7) and 30 646 of 30 646 with one 3-term relation (T, n = 6, 7), so 298 337 cases in all. That H-030's over-admission is exactly gate (no *path* in the kernel) vs injectivity (no *element* in the kernel) follows from Ladkani 2.3(c), together with the `isTiltingAtVertex` docstring and E-151 / E-318 (the gate does not imply J = 0). It does not follow from the sample, which only checks that two implementations of the cited test agree with a third (Hom). The test is not new: it is Ladkani 2.3(c), it is `tiltingPlus` (round 001), and the library has it as `mutation.isTiltingAtVertex` behind `search(..., refuseNonTilting=False)`. What is new is the shape split of the 8 and the agreement on this sample.
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

| sample (closed seed interval) | parents | gate F, Hom F | gate T, Hom F | gate T, Hom T | inj != Hom | shapes of gate T, Hom F |
|---|---|---|---|---|---|---|
| N5 [0, 5 999] | 319 | 770 | 1 | 1 593 | 0 | 1 core+1 |
| N5 [6 000, 35 999] | 1 673 | 3 996 | 8 | 8 314 | 0 | 6 exact, 2 core+1 |
| N6 [0, 3 999] | 617 | 1 678 | 7 | 3 779 | 0 | 5 exact, 1 core+1, 1 other (3277) |
| N6 [4 000, 39 999] | 5 167 | 14 153 | 34 | 31 300 | 0 | 18 exact, 10 core+1, 2 core+2, 4 other |
| N7 [0, 11 999] | 3 305 | 10 328 | 28 | 23 706 | 0 | 10 exact, 11 core+1, 7 other |
| P6 [0, 2 999] | 1 849 | 2 494 | 0 | 12 183 | 0 | - |
| *added in revision (`skeptic_inj2.py`, outputs `skeptic_inj2_<KIND>_<n>_<SEED0>_<COUNT>.out`)* | | | | | | |
| N8 [0, 9 999] (228 s) | 4 511 | 16 613 | 45 | 36 966 | 0 | 22 exact, 7 core+1, 2 core+2, 14 other (13 multi-arrow, 1 no input relation) |
| K5 [0, 35 999] (21 s) | 1 992 | 4 824 | 9 | 9 849 | 0 | 6 exact, 3 core+1 |
| K6 [0, 39 999] (98 s) | 5 784 | 16 078 | 40 | 34 833 | 0 | 23 exact, 11 core+1, 2 core+2, 4 other (3 multi-arrow, 1 no input relation) |
| K7 [0, 11 999] (98 s) | 3 305 | 10 477 | 24 | 23 561 | 0 | 9 exact, 10 core+1, 5 other (4 multi-arrow, 1 no input relation) |
| T6 [0, 39 999] (18 s) | 659 | 2 218 | 9 | 4 023 | 0 | 1 exact, 1 core+1, 7 other (1 multi-arrow, 6 no input relation) |
| T7 [0, 39 999] (99 s) | 2 222 | 8 604 | 29 | 15 763 | 0 | 1 exact, 1 core+1, 27 other (11 multi-arrow, 16 no input relation) |

In the added rows, gate F with Hom T is 0 everywhere. K has the same quivers and (vertex, side) set as N at the same n and seeds; only the binomial coefficients differ. T's "parents" counts only parents that received a 3-term relation. The classifier `e080` skips relations that are not binomial, so in T most kernels coming from the 3-term relation fall under "no input relation". "other" is broken down above by the reason the classifier prints, not by hand.

**How "other" is split (hand-read).** The classifier `e080` recognises only an *input* binomial relation of arrow length 3. It reports anything else as "other" with one of the reasons "k arrows on the side", "kernel at k sources" or "no input relation explains the kernel". It does not test consequence relations. Every finer split below is my hand reading of the printed `other:` lines plus the relations printed by `seeds`: "consequence" vs "two-arrow", which of two arrows kills and which identifies, and nn for 3277. The referee re-checked 7926, 23194, 696 and 3277. No command produces this split.

"other" at n = 6 beyond the 8: 23194, 33200, 37777 have two arrows on the side (23194: one arrow kills both p, q, the other identifies them); 7926 is one-arrow but the identifying relation is a consequence (1-2-3-6 = 1-4-5-6 = 1-2-5-6). At n = 7: 6 two-arrow cases (seeds 4774 6445 6740 9605 9813 11463) and 1 consequence case (696). "core+k" = E-080 mechanism with k extra surviving paths u ~> i. Gate F with Hom T: 0 in every sample (gate is necessary, as stated in `isMutable`).

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj.py seeds              # 1 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj.py planted            # 1 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj.py sweep N 6 0 4000   # 11 s (E-312 sample)
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj.py sweep N 6 4000 36000   # about 2 min
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj.py sweep N 7 0 12000  # about 2-3 min
```
(and `sweep N 5 0 6000`, `N 5 6000 30000`, `P 6 0 3000`, each under 2 min.) The last two arguments are SEED0 COUNT, so `sweep N 6 4000 36000` is seeds [4 000, 39 999].

Added in revision (`skeptic_inj2.py KIND n SEED0 COUNT`; KIND N = plain generator, K = non-unit binomials, T = one 3-term relation):
```
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj2.py N 8 0 10000    # 228 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj2.py K 5 0 36000    # 21 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj2.py K 6 0 40000    # 98 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj2.py K 7 0 12000    # 98 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj2.py T 6 0 40000    # 18 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj2.py T 7 0 40000    # 99 s
```
(times measured with 4 runs in parallel)

## Prior record

E-080 (the shape, one hand-built family), E-112 (lemma: for monomial + difference relations with scalar 1, J != 0 iff a circuit in Gamma_i; out-degree 1 gives p1 b = p2 b; nn and longer circuits not on walks), E-059 (2.3(c) vs gate on walk parents), E-126 (both-die, dim J = 1), E-312 (the 8 seeds, minimality), R-016. The injectivity test is `tiltingPlus` / `isTiltingAtVertex`, used in many entries (E-147, E-151, E-161: 1 081 step tests agree with an independent Hom test). So "an injectivity test closes them" is a recorded theorem applied to 8 more cases, not a discovery. New: the per-seed shape table (5 / 2 / 1), the observation that 3277 is an nn circuit (E-112: absent on walks; present among random parents), and 0 disagreements between 2.3(c) and Hom tilting on 298 337 random cases at n = 5-8: 114 372 N/P at n = 5-7, 53 624 N at n = 8, 99 695 with non-unit binomial coefficients, 30 646 with one 3-term relation.

## Code changed

None in the library. New scripts `workshop/rounds/126/skeptic_inj.py` and (revision) `workshop/rounds/126/skeptic_inj2.py` (imports the round-118/091 engines unchanged; signed coefficients as E-312). No tests touched.

## Next

- H-030 can be closed as *explained*: the gate is the path-wise shadow of 2.3(c); the fix exists (`refuseNonTilting`). The open question is only the library default (T10 request to me; cost E-151 8-11 %).
- theorist: E-112's lemma covers 3277 and 7926 only if "balanced" circuits through consequence relations are allowed; state it for relations given by a Groebner basis rather than as input.
- toolsmith: make `isMutable`'s docstring say that admission + `isTiltingAtVertex` is the tilting test, with seed 514 (exact) and 3277 (nn) as regression cases.
