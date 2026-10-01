# Experimentalist notebook (rewritten each round)

## What I now believe (after round 017)
- E-084 walk with the full-reduction `reduceAgainstPivots` (monkeypatch of 015 theorist_fix) gives identical tables
  for n = 7 c0-2, n = 8 c0, c1, n = 9 c0-3 (rejections 8/4/6, 4/14, 8/4/12/12; 0 guard-admitted failures).
  n = 8 c2: patched = unpatched through 16 500 expansions; the 10 M steps (key moved, gate+tilt) come after
  expansion 16 500 of 20 899 at depth 8, not reached by the patched walk in 10 min (17 058). Loose end unconfirmed by the
  walk; resting on E-085's replay. Files: `workshop/rounds/017/experimentalist_walk*.txt`.
- Patched run costs about 1.3x (Fractions). `MAXEXP=N` in `experimentalist_walk.py` gives a deterministic cap.
- From 014: n = 17 `5046`, `5056` orbits {0,2,4,6} = 122673, {1,3,5,7} = 54266 (closed). 4-letter words with a 4: big orbit = `444` orbit (sizes only compared).
- From 011: H-017 at n = 9 depth 6 negative (E-069); depth 7 overnight.

## What I tried
- 017: wrapper `experimentalist_walk.py` exec's `rounds/014/scholar_walk.py` after patching; runs in parallel 3-4 on 4 cores (inflates time);
  a time budget makes runs not comparable, so use MAXEXP for A/B.
- 014: `batch.py orbits 17` twice in parallel; orbscan4. 011: shards by `xargs -P 4`. 010: orbitclass n = 14..16.

## What I would do next
1. Overnight: n = 8 c2 depth 8 to completion (about 15 min single process, `--budget-sec 1500`), expect M = 0, rej 1; then depth 9.
2. Once the toolsmith's library patch lands, rerun the 017 rows with it (no monkeypatch), same files for comparison.
3. Row-set comparisons of n = 15 big orbit (3- vs 4-letter scan); 4-letter scan n = 16 by shards.
4. Depth 7 at n = 9 for H-017 (overnight), n = 17 key-coarser lists, 34x at n = 18, 44x (E-068).
- Watch: a timeout/budget is not a verdict; parallel timings are inflated; equal sizes are not equal sets; compare summary lines is weaker than comparing every rejection line.
