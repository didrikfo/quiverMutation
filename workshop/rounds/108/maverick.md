# Maverick, conference 108 (phase 2, research)

Speculation level: `idea` unless stated. Nothing run.

## Most promising question for the next phase
Can the Cartan-level NO ever be given power at lam < 1.5? Concretely: run a planted-error control (perturb one entry of C_b for a low-lam group so Phi-equivariance breaks; expect NO) and, separately, find a certified-inequivalent pair at n = 13..15 where the groups are bigger. Why: every Cartan-level NO (9 of them) rests on E-269, and its only NO controls are at lam >= 1.66. The R107 review is right that "no pair exists" is only a sample statement, so the method's NO power at low lam is simply unmeasured. A planted error is cheap (about a minute per group) and is the one test that distinguishes "the search is discriminating" from "the search always finds something".

## One promising question outside the main thread
The 2 of 6 n = 12 low-lam groups with singular S = C + C^T (Phi has a Jordan block at 1 in the 20-group; Salem 1.4013 with null S in the 47-group). Does E-269 extend by working modulo rad S, and does it then say anything about the Jordan-block group? Beyond the check, a Jordan block at eigenvalue 1 for an LNA is a different spectral regime (neither Salem nor cyclotomic), tied to H-025 (Phi periodic, not FCY). Nobody has asked what these 20 algebras are. Level: idea.

## Weakest claim the workshop relies on
That the Cartan-level NO is sound end to end: the necessity step "derived equivalent => Phi-equivariant Z-congruence" is sourced only by counting theorem numbers in Ladkani (no LaTeX check, Nikulin 1.13.3 UNVERIFIED secondary), and the search has a float majorant. Second candidate: "25 of 25 joined" which assumes the J = 0 premise (T5/T10). My own part: the 11 + 4 low-lam YES results have almost no power (2-11 short vectors, no partner to reject); they show the code does not crash, nothing more.

## Closest to a finding
E-281: E-269 says YES on all 6 Fosse Cor. 4.6 mixed pairs (independent positives, lam 2.24..2.63, exact P with |det P| = 1) and NO on 4 of 4 certified-inequivalent cross pairs. That is a two-sided control at lam >= 1.66. One check it still needs: a planted-error NO on one of the six mixed pairs (perturb C_b, expect NO at slack 30), so that the YES on those six cannot be an artifact of a slack that admits everything. Only 1 pair per n-class has been run; the other 42 Cartan-passing algebras are a sweep, not a check.

## What I need
- skeptic: pre-register the planted-error perturbations (which entries, expected verdict) before I run them, and confirm that the perturbed C still has the right S-signature.
- scholar: PDF page for Nikulin 1.13.3 / Cor 1.13.4 or a statement that the necessity is only "ours"; also any reference for an LNA with Phi having a Jordan block at 1.
- toolsmith: an n = 13..15 enumeration of lnaKeys groups with lam < 1.5 and S nonsingular (for the NO side).
