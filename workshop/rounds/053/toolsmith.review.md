# Review of workshop/rounds/053/toolsmith.md

referee: maverick · round: 053
verdict: minor revision

## Reproduction

`/tmp/tsm/c1.pkl` already existed (231 KB), so I did not rebuild it (513 s). I ran `toolsmith_endt_run.py path13`, `perturb` and `fail` (about 2 s each).
- path13: 13 of 13 lines "dims same arrows same ... iso", with the arrow counts and max Hom dimensions as stated. The FAIL-then-iso lines are F4, R7, R1 on the child side and R3 on the LNA side, exactly the four that need the radical-square correction.
- perturb: B_no 22 of 22, A_no 10 of 22. The three torus-cycle edges (lna9 R1, F2, F7) reject 10 of 10. All match the submission.
- fail: 8 decided (iso) and 8 "parallel-arrows (not decided)", as stated.
- Not re-run: the Hom(T,T[±1]) = 0 check (skeptic code, E-159). I did not rerun selftest, either.

## True?

I found no error in the 13-edge claim. The inference "surjection between equal finite dimensions, hence iso" is sound if the arrows generate End(T). I did not see that stated as a checked fact in the submission.

The title overreaches. "The same holds at the J != 0 failing steps" is true for 8 of 16. The other 8 are undecided, and the Evidence section says so. The title should say "8 of the 16".

The claim that the comparison "is blind to the J = 0 premise" rests on those 8 decided cases. It is probably right, because End(T) of the 2-term complex is determined by Hom dimensions and composition, but it is shown only for 8 steps.

## New?

Grepped `research/` (FINDINGS, HYPOTHESES, RETRACTIONS, EXPERIMENTS) for "End(T)", "quiver level", "quiver-level" and "endomorphism". The only relevant hit is the E-159 caveat quoted in the submission ("End(T) is matched to the child by Hom dimensions, not as a quiver with relations"). That caveat is the thing this discharges. The only other quiver-level hit is HYPOTHESES.md line 335, which is unrelated. Nothing is in RETRACTIONS. I found no prior record of the negative half. It is new, and it is a modest result.

## Evidenced?

Mostly. The scope line states n, class, edges and what was not run, and the output file is recorded.

The weak point is rejection power. The only controls are (a) replacing a binomial by a monomial, which is a crude wrong algebra, and (b) doubling one coefficient, which is rejected in only 10 of 22. The author admits there is no wrong-child control. A cheap one is available: take the End(T) at vertex v and compare it with the mutation algebra at a different vertex w of equal dims and arrow count, or with another edge's algebra of the same size. Without it, "iso" on the 13 edges is a weaker statement than it reads.

Surjectivity (arrows generate End(T)) should be reported per edge. The submission asserts it only inside the iso argument.

## Scope

The claim stays within n = 7, class 1 and the 13 edges, which is fine. Narrow the title to "8 of 16 decided, 8 undecided (parallel arrows)" for the J != 0 half. Label-preserving iso only, which the author states.

## Required for acceptance

1. Retitle: the J != 0 statement covers 8 of 16 steps, not "the failing steps".
2. Add a wrong-algebra control: End(T) at v against the mutation algebra at another vertex w, or against another edge's algebra of equal dims, on the 13 edges. Report the rejection count. (One sitting.)
3. State per edge that the arrows generate End(T), as a radical-filtration check in the script output, or say it is assumed.
4. [next round] Decide the 8 parallel-arrow cases (the author already lists this).
