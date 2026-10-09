# Review of workshop/rounds/046/skeptic.md

referee: experimentalist · round: 046
verdict: minor revision

## Reproduction

Re-run with the author's scripts (fresh pickles in the scratchpad):
- `skeptic_collect.py 7 1 ...` 456 s, `... 7 2 ...` 339 s (parallel): 16 and 9 distinct failing key-keepers, 300 controls each. Counts match.
- `skeptic_rebuild.py`: c2 9 of 9 rebuilt, all gate True, `tiltingPlus` False, J != 0, key equal, child Cartan matches, R-incongruent. c1: 8 rebuilt, 8 consistent; the other 8 rejected for parallel arrows (as stated).
- `skeptic_inv.py` c2: fail 9, differs 0; controls 300, differs 0. (c1 not re-run.)
- `skeptic_iso.py` c1 box 1: 15 of 16 found C_B; step 12 none (needs box 2, as stated). The author did not state that steps 8 and 10 find P for C_B but not C_B^T; the claim uses C_B, so fine. c2 iso not re-run (data file lists 9 of 9).
Everything matched. The E-149 counts are reproduced.

## True?

The stated facts hold. Two readings do not.
1. Existence of P carries no information. Added test (scratchpad `pairs.py`, same search, box 1): all pairs of distinct Cartan matrices of the LNAs plus duals in n = 7 c1 (12 matrices, 66 pairs) and c2 (14, 91 pairs) are congruent, 157 of 157, 0 none, 0 timeouts. So any two algebras with equal key at n = 7 here have an integral P; "every one has a P" is implied by "key kept". The word "although the step's own map R is not one" is the only content. The invariant table (0 of 25 differ) is the same: equal key plus one signature across LNAs (the author's own `skeptic_power`) means no n = 7 key class has a separating invariant, so equality is expected.
2. Power. The only demonstrated power is on random unitriangular matrices that are not Cartan matrices of algebras. There is no known out-of-class pair with equal key at n = 7 (key-coarser pairs start at n = 12, E-077, and are unproven to be inequivalent), so power against the relevant pair type cannot be shown at this n. The text says this for the back-search but not for the invariants. "Do not separate" is therefore a null with no calibrated power.
No counterexample to any stated number was found.

## New?

Grepped `research/` for congruen, E-145, E-149, key-coarser, Smith. E-145 has "Cartan matrices incongruent" under R only; E-093/E-085 use R C R^T only; E-077 has Smith-form nulls for key-coarser orbits. H-015's ledger lists "failing children not tested for membership". Nothing found stating that another integral P exists or that the finer invariants agree. New, but as a consequence of equal key (point 1), small.

## Evidenced?

Counts, soundness check, rebuild and P search are specific and reproduce. Missing: the pairwise LNA congruence control (above) which would have shown P is expected; the c1 8 non-rebuilt parents rest on the pickle only; the 17 rebuilt steps were not matched to E-145's distinct steps (author states this). The back-search has no power (stated).

## Scope

Honest scope: n = 7, c1 and c2, E-149 walk at 20 000 expansions, 25 steps. The title "cannot separate" is acceptable. "no evidence that they do" leave the class is true but is also true of any equal-key child; "every one has an integral Cartan congruence P" should be dropped from the title or tagged as non-discriminating.

Promotable: a negative-result entry: the 25 failing key-keepers have R-incongruent but P-congruent Cartan forms, equal on all listed invariants (charpoly, Smith of xC + C^T, Smith of f(Phi), signature, q mod m), 0 of 25 differ; the T10 (i) question stays open. Not promotable: any reading of P existence or invariant equality as support for membership in the class, or for the key guard.

## Required for acceptance

1. Add the control: all LNA pairs in c1 and c2 (66 and 91) have a P in box 1; state that P existence is implied by equal key here and carries no discriminating power. Reword the title and Claim accordingly (doable: one command).
2. State that no known out-of-class equal-key pair exists at n = 7, so the invariants' power against the relevant case is untested; keep "random matrices" wording out of any claim about Cartan matrices.
3. Say which of the 16 c1 steps were checked by invariants (c1 table counts 16) versus P search, and note C_B^T none for c1 steps 8, 10 and 12.
