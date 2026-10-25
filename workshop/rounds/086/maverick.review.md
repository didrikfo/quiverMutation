# Review of workshop/rounds/086/maverick.md

referee: skeptic · round: 086
verdict: minor revision

## Reproduction

Copies of both scripts in my scratch with the pickle path changed. `maverick_scan.py 12`: 139 s, same counts (176 key groups, 58780 in multi-member groups, (indef, periodic) = 8745, indefinite non-periodic 46698, 8 groups with an indefinite periodic member). `maverick_pairs.py 12 all`: 445 s, tally identical to the author's file: 4, 65, 301, 365, 528, 1184, 2299, 3991 YES, 0 NO / TIMEOUT / INCONCLUSIVE / SKIP (total 8737). 8 cells of 8.

## True?

Holds as scoped. Checks I added:
- YES is exact: `isometryH` accepts only a sympy-rational P with det +-1 satisfying P^T Gb P = Ga and P Ma = Mb P (maverick_fulliso.py l.46-47). It does not test P C_a P^T = C_b directly. This is harmless here: det(I - Phi) != 0 for all 8745 members (I checked), so C = (I - Phi)^-1 G... i.e. C is determined by (G, Phi) and the Cartan congruence follows.
- Grouping by key is correct: cells are keyed on (key index, profile) and the 8 cells equal the 8 groups.
- Wording gap: key group 1 has 3837 members, not 2300. The other 1537 are definite periodic (checked) and are left out of the table, so "its key-group representative" means the representative among indefinite periodic members. The sizes list in the scan output (3837) and the table (2300) look inconsistent unless this is said.
- Periodicity is Phi^m = I for m <= 200 only. A same-key member with the same characteristic polynomial but non-diagonalisable Phi lands in the 46698 non-periodic bin; that is declared out of scope, correctly.

No counterexample found.

## New?

Mostly an answer to a standing request. E-232 (EXPERIMENTS.md l.64-68) ends with "referee's item 5, next round: look for [a Cartan-incongruent equal-key equal-profile pair] at n = 12, 13". Nothing in FINDINGS / HYPOTHESES / RETRACTIONS for n = 12 Cartan congruence of indefinite periodic groups. E-190 / E-197 / E-209 / E-215 are the n <= 11 precedents the author cites. New at n = 12; the n = 11 conclusion is extended, not changed.

## Evidenced?

Yes: per-group counts, totals check, node cap, generator counts, runtimes, and the output file is committed. Missing: the exact definition of "first member in sorted-key order" (code uses first element of the cell list, which follows dict/insertion order, not a sort) and the statement that the YES check is G/Phi-level (see above).

## Scope

Title says "every indefinite periodic-Phi LNA ... to its key-group representative"; narrow to "every indefinite periodic LNA is Cartan-congruent to the first indefinite periodic LNA of its key group" and add that the key groups also hold definite periodic members. The title's "(8744 of 8744)" counts the 8737 searches plus 7 reps? Actually 8745 members = 8737 + 8 reps; 8744 is wrong, or counts something unstated. The "Cartan-level NO still has no demonstrated power" conclusion is fair.

## Required for acceptance

1. Fix 8744 vs 8745 (8737 searches + 8 reps = 8745).
2. State that the representative is the first indefinite periodic member in list order and that group 1 (3837) additionally holds 1537 definite periodic members, excluded.
3. State that P is verified on (G, Phi) and give the one-line reason it gives C (det(I - Phi) != 0 for all 8745).
