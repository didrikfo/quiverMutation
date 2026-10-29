# Theorist notebook (rewritten round 098)

## What I believe now
- 098 (T4): for the 24 floating width-5 rules, class A/B is decided by the first step alone via E-187/E-205 (24/24; B iff the first mutated arrow starts a relation of length >= 3, A iff length 2; lengths 6..8). Over all 414 rules at lengths w+1..w+2: 408 agree; 6 exceptions (w6: 2, w7: 4) leave at step 2. No B rule replays; E-187 is silent after the first step. Scripts: workshop/rounds/098/theorist_b5*.py, theorist_b_all.py.
- 095 (T4): H1's mechanism splits the 414 floating rules into class A (110; E-187/E-205 alone replay the rule: 0 disagreements at w+1..w+3, w+4 for w=6,7,9,10) and class B (304). For A, H1 at all lengths follows from a locality sketch (item 3 only by replay). Engine keeps labels: right at p swaps positions p,p+1; left p-1,p (empirical n<=8).
- 093 (T5): CONDITIONAL on omega-reduction of E-235: d_1 lands in Lemma A of E-240 (600 graded parents, 0 failures). Unproved: "d lowers ||.|| by 1", d^2 = 0 beyond E-202, Thm 7.5.
- 090: loops of M (Oppermann Steps 1-4): the 25 failing n=7 parents have acyclic Gabriel quiver, no loop at i.
- Older: repo step-7 kernel = Oppermann's iff D(i) acyclic (E-199); Gamma_red = N (E-217); LNA single-mutation rule proved (E-187, E-205).

## What I tried
- 098: wrote step-1 verdict vs engine; first run was wrong because applyAt returns the sequence shifted to the window start (d[3] is the unshifted rule). Lesson: take sequences from applyAt, never from the rule.
- 095: row-model replay of all 414 rules (labels, not positions). 093, 090, 089, 086, 083: see git history.
- Lessons: read the construction of the input before sampling it; check a model against the engine on multi-step sequences; suspicious uniform "gate" failures mean my harness, not the maths.

## Next
- T4: shape of the 14 width-5 middle algebras (r_v >= 3 case: Gabriel quiver, extra arrow?) and a one-step rule for them; if local, H1 is a theorem for B too. Prove the r_v >= 3 criterion beyond width 5.
- T4: the 6 step-2-fail rules (w6/w7) are the smallest test of a two-step rule.
- T4: class A width 11 at w+4 per rule; proof of label-swap rule from procedure step 7.
- T5: derive "d lowers ||.|| by 1" from Def. 3.2; d^2 = 0; Thm 7.5.
- Blind spots: proofs relative to printed v1; capped enumerations are not laws; r_v>=3 criterion observed at lengths 6..8 only.
