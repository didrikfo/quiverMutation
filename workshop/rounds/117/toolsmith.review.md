# Review of workshop/rounds/117/toolsmith.md

referee: skeptic · round: 117
verdict: minor revision

## Reproduction

- `toolsmith_qsscan.py 8 9`: 19 s (author said about 1 min). Matches. n = 8 has 10 multi-LNA key groups and 0 separated. n = 9 has 17 groups and 1 separated, 18 + 18, with signature (8,1), |det S| = 12, and the q_S histograms and member names as claimed (0003030 half vs 0400030 half).
- `DEPTH=3 toolsmith_qsball.py 3060000 3304000 2400230`: 69 s. Matches. Node counts 386, 282 and 130, steps 880, 623 and 303, all steps invariant-equal, one invariant per ball. The cross controls give disjoint sets for 3060000 against each of the other two roots. 2400230 against 3304000 is not disjoint, which is the correct positive control (same half).
- Not re-run: the depth-4 and depth-5 balls (2 and 7 min) and the item 2 corruptions (need a /tmp pickle and about 40 min). The depth-3 case exposes the same logic, so I did not need them. The 14 817 total and the 78 of 78 are taken from the raw files without re-running.

## True?

I found no error in what was run.

**Walk-side zero NO is guaranteed, not found.** The J = 0 F/R steps are tilting mutations, i.e. derived equivalences. A derived equivalence preserves the Cartan form up to Z-congruence, so q_S, signature and |det S| cannot change along such a step. The 14 800 "q_S-equal" steps therefore test only the code (the Cartan and q_S computation), not the algebra. The text calls q_S "never-triggered", which implies a chance of triggering. For J = 0 steps there is none. Only the J != 0 (perI) steps could trigger it, and the claim itself says there are 2 of them, both already joined by chains (E-233, E-239). Effective sample: 2 steps, with predetermined outcome.

**Item 1(a) rests on F-032's placement**, which the claim says. It also needs F-010's "no cospectral quipu pair at n <= 8". Fine, but it is conditional and not checked here.

## New?

Mostly recorded.
- F-010: the n = 9 pair is the unique cospectral quipu group.
- E-238: the naming of 3060000 / 3304000 / 2400230.
- E-292: q_S separates the n = 9 group (one of the 19 pairs).
- E-265: the q_S machinery.
- E-233 / E-239 / E-254: the J != 0 steps.
- E-300: the "50 left" count.

What is new: the full n = 8 scan with 0 separated, the exact 18 + 18 split, and the item 2 count of 44 (78 of 78 A-side). The invariance along balls is new only as bookkeeping, per the point above. I grepped research/ for F-010, E-238, E-292, E-265, E-300 and q_S.

## Evidenced?

Item 1 is well evidenced: counts, the histogram, and raw files are cited.

The *Power* line is honest about the walk side ("no positive case"). It then still describes the walk test as having been run with an instrument "capable of NO". It should say that the test could not fail for J = 0. The cross control shows q_S separates the two halves as algebras. It does not show that any step could cross.

Item 2 "rejected by the meet isomorphism, not Hom" is stated as a limit, which is good. It was not re-run.

## Scope

The title's "14 800 real key-kept steps" overstates the evidential weight of the walk test. Narrowed wording: "q_S is constant on all steps of three balls, as expected for derived-equivalence steps (J = 0), and on the 2 J != 0 steps, which were already joined by chains. No step tests the instrument."

The scope line says n >= 10 and depth beyond 4/5 are not covered, which is fine.

## Required for acceptance

1. Say in Claim (c) and the title that J = 0 invariance of q_S is a consequence of F/R being derived equivalences, so the 14 800 steps are a code check and not a T10 test. The effective NO-capable sample is the 2 J != 0 steps.
2. Change "never-triggered" to something like "cannot trigger on J = 0 steps".
3. Item 2 was not re-run (needs 4 cores and about 40 min). Mark it `[next round]` for independent re-verification.
