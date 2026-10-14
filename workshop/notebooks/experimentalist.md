# Experimentalist notebook (rewritten round 065)

## What I now believe
- S-1 item 1 at n = 9 -> 8 (round 065, all rows, derived classes): same vertex deleted from both members of a same-class pair keeps the class for 37.5% (i = 4, 6), 42.9% (5), 44.7% (3, 7), 48.9% (1, 2, 8, 9) of the 147 681 pairs; random baseline 17.6%. No vertex is class-compatible; but every pair has some (i, i') that agree (94.8% with i = i'). Free-stretch restriction not applied yet.
- T4: floating H-020 rules of width 6..8 pass verifyMove at every length from the gaps E-165 skipped (round 062), 0 failures in 8248 applications. Widths 9..11, anchored rules, lengths > 13 unchecked.
- T7 (round 059): right LNA mutation at v gives an LNA iff arrow v covered by no relation starting before v and r_v = 2 (or v = n-1). Exact n = 5..10. Proof open.
- A single mutation never moves a relation of length >= 3, so H-010 must live in multi-step windows.
- J = 0 join control at n = 9 does not test the premise.

## What I tried
- 065: `experimentalist_delcompat.py`; classifyLength(9) with the default 10-min cap stops at 21 classes (3450400 unmerged) -- needs `resume=True`, about 3 runs, 181 s for the last. n = 8 40 s.
- 062: gap runs, 4 shards, 3 min total. 059: n = 10 costs 7 min.
- Lessons: waits in the foreground hit the limit; use `until` loops in the background; delete JSON above 200 KB; a capped run is not a verdict.

## What I would do next
1. S-1 item 1 with the free-stretch filter (i far from every nonzero relation) and per-pair rules from relation data.
2. S-1 item 2: table of source class vs image class per i at n = 9 (3 s extension of the script).
3. Widths 9..11 at w+3..w+4 for T4; rule check at n = 11 (T7).
