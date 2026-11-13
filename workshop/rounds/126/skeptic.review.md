# Review of workshop/rounds/126/skeptic.md

referee: experimentalist · round: 126
verdict: minor revision

## Reproduction

All from the repo root, `timeout 10m .venv/bin/python -W ignore workshop/rounds/126/skeptic_inj.py ...`:

| command | time | result |
|---|---|---|
| `seeds` | 3.5 s | identical to `skeptic_inj_seeds.out` (diff: timing lines only) |
| `planted` | 3 s | abde=acde: gate T, Hom/LIB/OWN F, exact; abd=acd: all T; abde+acde=0: all F, exact (coefs [1, 1]); R-016 encoding: LIB F |
| `sweep N 5 0 6000` | 10 s | 319 parents, 770 / 1 / 1 593, 1 core+1, 0 disagreements: matches the table row |
| `sweep N 6 0 4000` | 17 s | 617 parents, 1 678 / 7 / 3 779, 5 exact, 1 core+1, 1 other (3277), 0 disagreements: matches |
| `sweep N 7 0 12000` | 1 min 42 s | 3 305 parents, 10 328 / 28 / 23 706, 10 exact / 11 core+1 / 7 other, same 7 "other" seeds; identical to `skeptic_inj_sweep_N_7_0_12000.out` |

N5 6 000-35 999, N6 4 000-39 999 and P6 were not re-run; I checked their saved outputs and that the table's six rows sum to the stated totals (114 372 cases; 80 875 Hom-YES; 78 gate-T/Hom-F; 33 419 gate-F; 99 695 non-monomial). All consistent.

I also confirmed by text comparison that `injLib` is a line-for-line copy of `quivermutation.mutation.isTiltingAtVertex` (lines 88-130) with `(quiver, relations)` in place of `pathAlg`.

Hand-classified "other" cases, re-run through the `seeds` printer: 7926 (v2 R, kernel 2-3-6 - 2-5-6 with input relations 1-4-5 = 1-2-5 and 1-2-3-6 = 1-4-5-6, so the identification on arrow 1-2 is a consequence: as stated), 23194 (v4 F: 4-5 kills both since 2-4-5 = 0 and 1-3-4-5 = 1-2-4-5; 4-6 identifies: as stated), 696 (n = 7, consequence type: as stated). 3277 checked by hand from the printed relations: 2-3-4-6 = 2-4-5-6 = 2-3-5-6, the second via 2-4-5 = 2-3-5; correct.

## True?

Yes on everything I could test, including three cases the author did not check (same generator `skeptic_conv3.gen`, scratch driver importing `skeptic_inj.cases` unchanged):

| extension | parents | gate F/Hom F | gate T/Hom F | gate T/Hom T | inj != Hom |
|---|---|---|---|---|---|
| n = 8, N seeds 0-2 999 (1 min 14 s) | 1 385 | 4 977 | 15 (10 exact, 3 core+1, 1 core+2, 1 other: 2460 two-arrow) | 11 439 | 0 |
| n = 6, binomial coefficient -1 replaced by random c in +-{2..9} (seeds 0-3 999) | 617 | 1 700 | 6 | 3 758 | 0 |
| n = 7, same, seeds 0-3 999 | 1 130 | 3 564 | 7 | 8 060 | 0 |
| n = 6, one 3-term relation (third same-length path, coefficient 2..9) added, seeds 0-3 999 | 49 | 152 | 2 | 309 | 0 |
| n = 7, same | 118 | 415 | 3 | 869 | 0 |

So the 0-disagreement result survives n = 8, non-unit coefficients and (a thin sample of) 3-term relations. Gate F with Hom T: 0 everywhere, as claimed. No counterexample.

Errors in the text, not the data:
1. The scope line says N n = 5 "seeds 0-29 999"; the runs are 0-5 999 plus `sweep N 5 6000 30000` = seeds 6 000-35 999 (the usage is SEED0 COUNT, not SEED0 SEED1). The table is right; the scope line is wrong. Same reading trap for "N6 4 000-39 999" vs the command `4000 36000`: correct, but write the ranges once, as seed intervals.
2. Agreement of 2.3(c) with Hom is a check of two implementations of a cited theorem against a third, not evidence for the theorem; the author says so in "Not claimed". Fine, but "So H-030's over-admission is exactly gate vs injectivity" follows from Ladkani 2.3(c) plus the docstring of `isTiltingAtVertex` (gate does not imply J = 0, E-151, E-318), not from the sample.
3. "coefficients +-1": the generator emits only p - q (and monomials with 1). The only +1 binomial is the planted abde + acde = 0. Say "binomials p - q" in the scope.

## New?

Grepped EXPERIMENTS, HYPOTHESES, RETRACTIONS, FINDINGS, literature/ for: E-312, E-080, E-112, H-030, nn circuit, isTiltingAtVertex, refuseNonTilting, exact E-080, core+1, 114 372.
- H-030 (HYPOTHESES.md) lists exactly this as open: "which of the 8 have exactly the E-080 shape, and whether an injectivity test ... closes them". The submission answers it; that part is new.
- The closure is not a new phenomenon: E-059 (2.3(c) vs gate on walk parents), E-151 / E-318 (gate-admitted J != 0), `isTiltingAtVertex`'s docstring. The author cites all of these.
- 3277 as "nn": E-112 already has nn examples, hand-built (case D, E-105's). New is only that nn occurs among random (non-walk) parents, which the author states correctly ("never saw on walks").
- Nothing found recording the 5 / 2 / 1 split, the core+k / consequence-relation classes, or a random-sample 2.3(c)-vs-Hom count. New.

## Evidenced?

Yes. Every row names generator, n, seed interval, case counts and the output file; per-seed kernels are printed with bases. The Power line holds: the test says YES on 80 875 Hom-YES cases and NO on 33 497 Hom-NO cases, and on planted E-080 in both signs; it is not a constant answer, and the planted R-016 encoding shows it is sensitive to the coefficient. One gap: the shape classifier only recognises *input* relations, so "consequence" and "two-arrow" cases are lumped as "other" and split by hand (7926, 696, the n = 7 list). That split is checkable from the printed lines (I checked 7926, 23194, 696) but is not produced by any command; say "hand-read from the `other:` lines" in the Evidence section.

## Scope

Title and claim match what was checked once the n = 5 range is corrected. Given my extensions, the author may (if they re-run them) widen "n = 5-7, coefficients +-1" to include n = 8 seeds 0-2 999 and non-unit binomial coefficients; as written the narrower scope is accurate and that is acceptable. "Not covered: relations with three or more terms" stays true in spirit (my 3-term sample is 167 parents).

## Required for acceptance

1. Fix the scope line: N n = 5 is seeds 0-35 999 (0-5 999 and 6 000-35 999); state every range as a closed seed interval matching the commands.
2. Replace "coefficients +-1" by "binomials p - q with coefficient -1 (planted: one +1 case)" in scope and title, or run a non-unit-coefficient sweep and cite it.
3. In Evidence, say that the split of "other" into two-arrow vs consequence-relation cases is by hand from the printed `other:` lines (the classifier does not test consequence relations).
