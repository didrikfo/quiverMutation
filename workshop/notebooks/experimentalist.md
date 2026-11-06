# Experimentalist notebook (rewritten round 119)

## What I now believe
- T7 / H-010: no gap-2 pair has ever lowered overlap (E-288, round 119); runs of three do (2->0, 3->2). The two-bystander shape has no control. k >= 5 is hours per cell: 109 s for (8,4)(10,4) at k = 4 margin 1, capped at k = 5, 6. Undecided is not negative.
- Harness: at k = 5, 6 `arrowPaths.liftPath` raises "parallel arrows does not name a path" for (8,3)(10,3). Deeper cells may need that fixed.
- S-1/H-026 (round 117): at even n the law K0 = min(h, n-4-h) holds off-centre; the self-mirror class has no split (n = 10, 12, 14; 16 K >= 3). Needs rewording; refinement fitted post hoc. Rivals min(h, n-5-h), min(h, n-3-h) are equal off-centre. Class level (n = 15 h = 3): image key constant on forward-rewrite components only.
- Two-3 classes: no law (115). E-263 (4,7) was a writer bug (R-015).
- Practical: kpower.py `table N K1 h,h,..`; `orbits` uses lnaMoves.rewritesOf. Lesson: if every fit is at one parity, test the rival that differs at the other.

## What I tried
- 119: `rounds/119/experimentalist_gap2.py` (+ _k4/_k5/_k6 outputs). 117 kpower; 115 kloc, nonlone; 114 tail; 113 kthresh; 110 refusal; 107 chain.

## What I would do next
1. T7: overnight k = 5, 6 gap-2 pairs (see submission 119), or ask theorist for an argument that gap-2 cannot lower.
2. Self-mirror exception at even n; non-lone self-mirror classes.
3. Class level with q_S / Z[Phi] at n = 15 h = 3.
