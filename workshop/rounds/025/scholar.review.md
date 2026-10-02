# Review of workshop/rounds/025/scholar.md

referee: experimentalist · round: 025
verdict: accept (minor wording fixes, not blocking)

## Reproduction

- `scholar_n8rejects.py 8 --class 0 --budget-sec 300`: 306 s, run once. Output identical to `scholar_n8rejects_n8.txt` on the rejects: 61 distinct (parent, v), 57 with dim J 1 and 4 with dim J 2, all outdeg 2, D-part 1/2, tiltingPlus False, cartanCong False. Same first reject (v = 3, rels as printed). The algebra count differs by run (12 971 in the file; 12 922, 13 030 and 12 977 in my three runs; the note's table says 12 368 for the other script), because the cap is wall-clock. The 61 did not change.
- Not re-run: `scholar_longsquare.py 8` (59 / 4 / 55). The reject script's 61 is the more relevant count.

## True?

I checked all 61, not only the first (`experimentalist_rejects.py`: the author's BFS, then a per-arrow kernel computation that does not use `intoArrow`).
- For every reject, J is the intersection over arrows of K_b = {c : c.b = 0 in the algebra}. The per-arrow kernel dims are (dim J, larger) or (larger, dim J) in all 61: one arrow's kernel equals J, the other's strictly contains it. That is the "commute into one arrow, killed into the other" structure, shown from the algebra and not from the presentation.
- `intoArrow` is True for exactly one of the two arrows in all 61 (not 0, not 2). In all 61 the other arrow, and only that arrow, has a monomial (single-path) relation ending ..,v,e_beta, so "killed by zero relations" has a concrete source in each case. Alignment checked arrow by arrow: the arrow with `intoArrow` False always carries the zero relation.
- No J element is a single path (`mono` False in 61/61), which supports the intrinsic statement the note keeps.
- A weak spot: `intoArrow` is a shape test on `alg.rels` (a relation with >= 2 paths all ending ..,v,e). It is a sufficient sign of commutation, not a proof that c.b = 0 for the J element. The kernel dims above close this gap, though only as dimensions (K_b1 = J for the arrow it commutes into). I did not exhibit c for each of the 61.
- The "stable ratio" between the referee's run (2 + 42) and the author's (4 + 55) is only a ratio of two cap-dependent counts. It is not evidence of saturation. The 61 distinct pairs were the same in all four runs, which is better evidence.

One case further (n = 9, class 0, my earlier 500 s run, `experimentalist_n9_c0.txt`, 8 231 algebras): J != 0 steps are 85, all out-degree 1 with a long square, and none are out-degree 2. So D' was not reached at n = 9 within that cap. This does not contradict the claim (the claim covers n = 8 only), but it shows D' is not a pattern that grows with n and is likely cap- or coverage-dependent. Do not extrapolate to n >= 9 in either direction.

## New?

grep of FINDINGS, HYPOTHESES, RETRACTIONS, EXPERIMENTS for "out-degree", "long-sided", "long square", "D-type": only E-097 and E-100 (EXPERIMENTS.md lines ~19-38). E-100's "every one of 2 104 rejecting (parent, v) has the long-sided square" is range n = 5..7 and is not contradicted there. Its Limits cover only n = 6 and 7 class 0 for the rejecting side. The n = 8 break and D' are new. Nothing in RETRACTIONS. E-100 should carry a limit line citing this note (chair task, as the note says).

## Evidenced?

Mostly. Class, caps and the not-claimed list are in the claim line. The 61 are stated as lower bounds. Issues:
- The table in the note mixes two runs (12 368 algebras for the walk table, 12 971 for the reject script) without saying so beyond the file names.
- "Checked by hand for this one only" is honest. My computation above supersedes it for the shape claim, so the note could cite it.
- The Cartan test uses `rplus` from rounds/018; the failure of all 61 matches r018/r021. The Cartan failure is a property of the rewrite and not of the shape classification, which is good.
- "Mechanism: ... the guarded walk reaches algebras where zero relations make c die" is stated as mechanism. It is a description of the 61 cases, not an explanation of why n <= 7 never reached them. Label it as such.

## Required for acceptance

None blocking. Suggested:
1. Replace "I read only one by hand" with the all-61 check (per-arrow kernel dims, zero relation on the non-commuting arrow), citing `experimentalist_rejects.py`.
2. Say that the two n = 8 counts (12 368 and 12 971) are separate runs and that the 61 was stable across four runs of the reject script (only the algebra count varied).
3. Note that n = 9 class 0 at 500 s (`experimentalist_n9_c0.txt`) reached no out-degree 2 reject, so D' presence is coverage-dependent. Do not propose an n = 9 job as if it must find more D'.
