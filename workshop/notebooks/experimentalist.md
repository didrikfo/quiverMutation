# Experimentalist notebook (rewritten round 062)

## What I now believe
- T4: floating H-020 rules of width 6..8 pass verifyMove at every length from the gaps E-165 skipped (w6 L11, w7 L10-11, w8 L11-12): 0 failures in 8248 applications (round 062). Remaining unchecked: widths 9..11 beyond w+2, anchored rules at gap lengths, lengths > 13.
- `theorist_rulelen.py` has an `L >= w+5` filter that excludes gap lengths; use `workshop/rounds/062/experimentalist_gaps.py L W [START STEP] [--plan]`.
- T7 (round 059): right LNA mutation at v gives an LNA iff arrow v covered by no relation starting before v and r_v = 2 (or v = n-1); image = slide of the length-2 relation. Exact n = 5..10. Proof open.
- A single mutation never moves a relation of length >= 3, so H-010 must live in multi-step windows.
- J = 0 join control at n = 9 does not test the premise; do not quote as support.

## What I tried
- 062: gap runs, 4 shards each, total about 3 min on 4 cores (fast; I overestimated cost).
- 059: `experimentalist_which.py`, `_rule.py`, `_capcost`-style runs; n = 10 costs 7 min.
- Lessons: loops with `wait` run into the 120 s foreground limit, run in background and poll; delete JSON above 200 KB.

## What I would do next
1. Widths 9..11 at w+3..w+4 (54 rules), same script; likely cheap.
2. Rule check at n = 11 (T7), sharded.
3. Anchored rules at the gap lengths.
