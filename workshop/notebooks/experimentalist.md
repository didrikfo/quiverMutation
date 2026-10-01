# Experimentalist notebook (rewritten each round)

## What I now believe (after round 019)
- T5, Cartan vs tiltingPlus (E-093 follow-up, 019): on guarded walks n = 5..7 (n = 6, 7 class 0 only, 480 s caps) and the E-078 family, dim ker is 0 on every
  tilting step (0 of ~167 000) and >= 1 on every non-tilting step (0 of 1 056 with dim ker 0). Non-tilting totals: 1 on all 1 050 walk steps
  (n = 6: 907, n = 7: 143, each a distinct parent with exactly one bad vertex); 2 and 3 only in E-078 (3 / 2 / 1 steps at total 1 / 2 / 3). Per-i dim ker never > 1.
  Files: `workshop/rounds/019/experimentalist_kerhist*.{py,txt}`.
- From 017: E-084 walk with the fixed `reduceAgainstPivots` gives identical tables (n = 7, 8 c0-1, 9); n = 8 c2: loose end rests on E-085's replay. `MAXEXP=N` in `rounds/017/experimentalist_walk.py` is a deterministic cap.
- From 014: n = 17 `5046`/`5056` orbits {0,2,4,6} = 122673, {1,3,5,7} = 54266 (closed). 4-letter words with a 4: big orbit = `444` orbit (sizes only).
- From 011: H-017 at n = 9 depth 6 negative (E-069); depth 7 overnight.

## What I tried
- 019: copied scholar's script, added kerdims on every step and a parent tally; four runs in parallel on 4 cores (n = 6, 7 each 8 min wall cap).
  Time caps make step counts differ between runs of the same script (696 vs 907): quote them as lower bounds, not constants.
- 017: wrapper exec'ing `rounds/014/scholar_walk.py` after a monkeypatch. 014: `batch.py orbits 17`. 011: shards by `xargs -P 4`.

## What I would do next
1. Classify the 1 050 rejecting parents by shape (A5 or not): is dim ker <= 1 just A5? Then n = 6 classes 1-3, n = 7 classes 1+ (use a step cap, not wall clock).
2. Overnight: n = 8 c2 depth 8 to completion (`--budget-sec 1500`), expect M = 0; then depth 9.
3. Row-set comparison of n = 15 big orbit (3- vs 4-letter scan); 4-letter scan n = 16 by shards.
4. Depth 7 at n = 9 for H-017 (overnight); n = 17 key-coarser lists; 34x at n = 18, 44x.
- Watch: a timeout/budget is not a verdict; parallel timings are inflated; equal sizes are not equal sets; a hand-built family (E-078) is not the walk's distribution.
