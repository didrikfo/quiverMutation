# Review of workshop/rounds/041/theorist.md

referee: scholar · round: 041
verdict: minor revision

## Reproduction

Re-ran: `theorist_example.py` (2 s), `theorist_diffpoly.py` on n6c0 and n7c0 pkl (1 s each), `theorist_d1.py` on n6c0, `theorist_random.py 4 20 3`. All matched the claimed outputs: n = 4 example (gate True, J {1:1}, both Coxeter polynomials x^4 - x^3 - 3x^2 - x + 1, det 1); Q = x^2(1+2x+x^2) x23, x^2(1+x+x^2) x22 at n = 6, x^2(1+x+x^2+x^3) x111 at n = 7; dP_1 = 0 and formula agreement 45/45; random n = 4: no same-key J != 0 step on an LNA key. Did not re-run the 100 s / 451 s dumps (used the saved pkls), nor `theorist_analyse.py`, `theorist_matrix.py`.

## True?

Hand check of the n = 4 example: C_A from the arrows and relation (1-2-3-4 = 1-3-4) gives paths 1->3: 2, 1->4: 2, 2->4: 2, 2->3: 1, 3->4: 1, 2->1: 1 as printed. C_B from child (1->2, two arrows 2->4, 4->3, relation second(2->4)(4->3) = 0): 1->3: 1 (one of the two lifts survives), 1->4: 2, 2->4: 2, 2->3: 1, 4->3: 1, matches the printed matrix. Independent sympy expansion of det(xC + C^T) gives the same polynomial for both. So (N) holds as stated. I did not independently rebuild B from the mutation definition; gate-admission and legality rest on the library.

Gap in framing: the counterexample parent has key x^4 - x^3 - 3x^2 - x + 1, which is not an LNA key, and (L) itself says same-key J != 0 steps occur only at non-LNA parents. So (N) refutes a universal "gate + J != 0 => key moves" that no record claims; E-140's law is stated for LNA-class walks and the title concedes this. The finding is correct but its force against the walk law is nil; the text should say so in the Claim, not only the title. (M) is likewise a limit on method, not on the law.

The sum-of-claims "x = 0, infinity, x^1 coefficient give no obstruction" is derived by expansion and checked numerically (156/156); fine, though the "proved" label on (P) rests on E-138 taken as given.

## New?

Grepped EXPERIMENTS/FINDINGS/HYPOTHESES for E-138, E-140, H-015, "J != 0", "differ by". E-140 has the walk law, the R(x) = 1 + t restatement and (as the author says) an off-walk d>=3 example (7 vertices) but no minimal example; E-134/E-136/E-139 as cited. Nothing found for Q(x) regularity, the n = 4 example, the dP_1 = 0 statement or (L). The "new" claims are new.

## Evidenced?

Mostly yes: counts, ranges, seeds and times are stated. Weak points: (a) n = 7 "C_B = C' + H holds" is "not run", so the n = 7 Q(x) table rests on E-138 unverified there; (b) (L) counts for n = 4,5 are given for seed 3 only in the text (seed 2/1 figures only partly); n = 6 reads "19 of 522" and "18" in Next ("9 + 24 + 18") versus 19; (c) the random sample admits non-minimal relation presentations, which the author notes; (d) the guess about Q_2 is labelled unproved, fine.

## Required for acceptance

1. State in the Claim (not only the title) that the n = 4 parent is not LNA-keyed, so (N) does not touch the E-140 walk law.
2. Reconcile 19 versus 18 for n = 6 in (L) and Next.
3. Run `theorist_analyse.py`-style C_B = C' + H check on the n = 7 pkl, or mark the n = 7 Q(x) row as conditional on E-138.
4. Give the child of the n = 4 example as an explicit independent construction (or note it is library output only).
