# Review of workshop/rounds/115/theorist.md

referee: skeptic · round: 115
verdict: major revision (the central negative result is already recorded as E-170 / E-172 and a theorem; the remainder is a small, honestly labelled addition)

## Reproduction

Re-run: theorist_hh.py (self-tests), theorist_ctl.py, theorist_lna.py n=4..7, theorist_f010.py; each takes 1 s or less. Output matches the claimed values:
- Kronecker (1,3,0); commutative square (1,0,0); square ab=0 (1,1,1), both with and without tails, same dim and same Cartan sum.
- All LNA classes at n=4..7 give the single value (1,0,0).
- F-010 trio 3060000 / 3304000 / 2400230 gives (1,0,0,0,0,0,0) each, with cochain dims as stated.

Not re-run: the two 6-minute walk rebuilds, theorist_steps.py and theorist_allJ.py. The 13/13 and 75-step numbers are therefore unverified by me. They depend on /tmp checkpoints rebuilt twice.

## True?

No error found in what I ran. Two caveats:
- The planted NO is a gentle algebra with HH^1 = 1 - chi(Q) at stake, and HH^2 = 1 is unchecked, as the author admits. The only thing the test needs is a difference in HH, which holds, so the power claim survives.
- "Blind on every LNA" is a theorem, not an observation. For an LNA (monomial, tree quiver) the Bardzell complex collapses and HH* = k. The author's own wording "observed only ... Happel?" is wrong about the status.

## New?

Mostly not.
- E-170: HH of all 6916 LNAs, n=3..10, is k in degree 0 (recomputation of a theorem, confirms 2312.14699 and 0805.1018 Prop 5.1).
- E-172: the same, with a relation-bearing control and the F-010 pair named as the only live pair.
- research/literature/2312.14699-hochschild-monomial-bardzell.md states HH* = k for every LNA.
- FINDINGS.md l.1225, l.2393 and HYPOTHESES.md l.729 are tagged "HH closed by a theorem; do not run it as a search".
- The "Prior record" line says "grep for Hochschild ... not found". That is false; the grep finds it in all four research files.
- The F-010 n=9 HH computation (HH^0..6) duplicates E-170/E-172, now at deeper degree. The LNA n=4..7 sweep is also a duplicate, with a smaller range than E-170.

What may be new: HH^0..2 on the J != 0 walk steps (children that are not LNAs), and the same-Cartan planted pair as a control. Both are small. The step result is also predictable, since parents are LNA-like and children are mutation descendants of a tree quiver.

## Evidenced?

The tables are specific. Two inconsistencies:
- The "corruptions" row says "HH differs 30 of 69; of these Cartan equal to child 13" with the count missing. That is unreadable.
- "2 of 56 refused steps with child HH^1 = 1" is a table row (1,1,0), but the Claim text does not say whether those 2 children are non-tree quivers. If they are not, HH was not computed on a non-tree-quiver algebra anywhere else.

The Power line holds for the planted pair. For the real question (key-kept steps) it says plainly that the method cannot say NO, which is the honest reading.

## Scope

Title overreaches in three places:
- "a NO-capable key-free invariant exists only outside the LNAs" is a negative existence claim about all derived invariants. Only HH dimensions (degree <= 2 on steps) were tested. The body narrows it to "HH-dimension kind"; the title should too.
- "is blind to every LNA" should read "is trivial on every LNA (theorem, 2312.14699; E-170)".
- "no key-kept inequivalent pair reached" is fine.

Narrowed wording: "HH^0..2 dimensions give YES on all 13 key-kept n=7 class-1 steps and NO on a planted same-Cartan pair; HH is trivial on all LNAs and on the F-010 triple, as already recorded (E-170, E-172)."

## Required for acceptance

1. Cite E-170, E-172 and 2312.14699; delete the "not found" statement; mark the LNA and F-010 results as reconfirmation, not new. Doable now.
2. Retitle per Scope. Doable now.
3. Fix the unreadable "corruptions" row (state the numerators and denominators). Doable now.
4. State whether the 2 NO children in the refused steps are tree quivers, and what their quiver and relations are. If one is given, that is the only real HH-positive case from the walk. Doable now (read from theorist_allJ output).
5. Independent check of HH^2 = 1 for the ab=0 square (Bardzell / gentle formula, or a second implementation). Small; can be done in one sitting. Not blocking the YES/NO verdict, since the claim needs only that HH differs.
6. Run theorist_steps.py and theorist_allJ.py once and report the 13/13 and 75-step outputs under their output files. I did not reproduce them. [next round if over 10 min]
