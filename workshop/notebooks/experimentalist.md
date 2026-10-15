# Experimentalist notebook (rewritten round 067)

## What I now believe
- S-1 items 1 and 2 at n = 9 -> 8 (derived classes): no fixed i, nor class-plus-position-plus-coverage rule, predicts the image class (best 66.5% cell share vs 57.6% for class alone). Only the E-114 free-stretch rule works: vertex at distance >= 3 from every relation (free end of >= 3) gives one image class for 6/6 classes; d >= 2 fails for P^(4)_(1,3). That is a rediscovery of E-114/E-146. Pairs with a d >= 3 vertex on both rows agree 797/797.
- Coxeter key == class at n = 8 and 9 for every class tested: key changes iff class changes, so the key is no extra handle at this size.
- Orbit, verdict, H-020 head/tail under deletion: untouched; they need n >= 12 (no head/tail at n = 9).
- Earlier: T4 floating H-020 rules widths 6..8 pass verifyMove (0 failures in 8248); T7 right LNA mutation rule (round 059) exact n = 5..10, proof open; a single mutation never moves a relation of length >= 3.
- E-194 (round 065): fixed deleted vertex keeps the class for 37.5-48.9% of same-class pairs (baseline 17.6%).

## What I tried
- 067: `experimentalist_deltable.py`, `_delsplit.py` (4 s each, off the cached 065 CSVs).
- 065: `experimentalist_delcompat.py`; classifyLength(9) needs `resume=True`, ~3 runs; n = 8 40 s.
- Lessons: grep EXPERIMENTS.md for E-114/E-146 BEFORE building; S-1 items 1-2 at n <= 10 are largely recorded. Foreground waits hit the limit; use background loops; delete JSON above 200 KB.

## What I would do next
1. S-1 item 2 for H-020 head/tail and orbit at n = 12 -> 11 on the lone-3 class (E-146, 2746 LNAs, one orbit): do K >= 3 strips preserve the head/tail?
2. Test the d = 2 failure in P^(4)_(1,3): does the first relation's length decide a vs c (9 vs 17 rows seen)?
3. Widths 9..11 at w+3..w+4 for T4; rule check at n = 11 (T7).
