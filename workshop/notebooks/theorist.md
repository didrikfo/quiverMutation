# Theorist notebook (rewritten round 095)

## What I believe now
- 095 (T4): H1's mechanism splits the 414 floating rules into class A (110 rules; every step gated and an LNA, so E-187/E-205 alone replay the rule: 0 disagreements at lengths w+1..w+3, and up to 13/11/11/12 for w<=5/6/7/8) and class B (304 rules; first step leaves the LNA class, E-187 silent). For A, H1 at all lengths follows from locality: the replay reads only the window pattern, r_{s+w}, and end distances, so lengths w+1..w+4 suffice (item 3, image at r_{s-1}, is covered only by replay). Engine keeps labels: right at p swaps positions p,p+1; left p-1,p (empirical, n<=8). Row model = engine on 9 396 cases.
- 093 (T5): CONDITIONAL on omega-reduction of E-235: d_1 lands in Lemma A of E-240 (600 graded parents, 0 failures). Remaining unproved inputs of E-235: "d lowers ||.|| by 1", d^2 = 0 beyond E-202, Thm 7.5. Only failure mode: a loop at i.
- 093 (T10): Hom(T,T[-1]) at S=-1 is sum J_i (E-128 L2), a restatement of J != 0; a detector must name S=-1 and cannot be a NO.
- 090: loops of M (Oppermann Steps 1-4): four sources; the 25 failing n=7 parents have acyclic Gabriel quiver, no loop at i.
- Older: repo step-7 kernel = Oppermann's iff D(i) acyclic (E-199); Gamma_red = N (E-217); LNA single-mutation rule proved (E-187, E-205).

## What I tried
- 095 response: class A now also checked at w+4 for w=6,7,9,10 (0 failures, no rule changes A/B). Title now 'supported by replay'; label-swap rule flagged empirical in the Claim.
- 095: row-model replay of all 414 floating rules; first run gave 392 spurious gate failures until I noticed the engine relabels vertices after a mutation. Lesson: sequences name labels, not positions.
- 093: step-4 hypothesis check; Hom table. 090, 089, 086, 083 and earlier (see git history of this file).
- Lessons: read the construction of the input before sampling it. grep for the identity that makes a "detector" tautological. Check that an abstract model reproduces the engine on multi-step sequences, not only single steps.

## Next
- T4: one-step rule for LNA-with-one-extra-arrow (middle states of class B); smallest: 14 width-5 B rules. If local, H1 is a theorem for all rules.
- T4: class A width 11 at w+4 (w=9, 10 done in R095 response, split A/B stable at w+4 for w=6,7,9,10; w=11 exceeds 10 min, run per rule); proof of the label-swap rule from procedure step 7.
- T5: derive "d lowers ||.|| by 1" from Def. 3.2; d^2 = 0; Thm 7.5; apply Lemma A to gr_G of a real kM (n <= 5).
- Blind spots: proofs relative to printed v1; capped enumerations are not laws; the locality argument is a sketch, items 1-3 not machine-proved beyond replay.
