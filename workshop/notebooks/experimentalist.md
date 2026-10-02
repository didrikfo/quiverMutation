# Experimentalist notebook (rewritten each round)

## What I now believe (after round 022)
- T5 control (022): at n = 5..7 guarded walks (n = 6, 7 class 0, n = 6 class 1, n = 5 classes 0, 1) the long-sided square of E-097 occurs at 0 of 479 761
  tilting steps and at 100 % of rejecting (parent,v) (1 842 at n = 6, 262 at n = 7). Strict A5 is not selective: 0.5-1.3 % of tilting steps have it
  (807 / 150 592 at n = 6). So long square <=> rejection on the walks (conjecture, loose test, capped runs). Files: `workshop/rounds/022/experimentalist_shapectl*.{py,txt}`.
- T5 (019): dim ker is 0 on every tilting step and >= 1 on every non-tilting step; totals 1 on all 1 050 walk steps; 2 and 3 only in E-078. (`rounds/019/experimentalist_kerhist.py`.)
- From 017: E-084 walk with fixed `reduceAgainstPivots` gives identical tables (n = 7, 8 c0-1, 9); n = 8 c2 loose end rests on E-085's replay. `MAXEXP=N` in `rounds/017/experimentalist_walk.py` caps deterministically.
- From 014: n = 17 `5046`/`5056` orbits {0,2,4,6} = 122673, {1,3,5,7} = 54266 (closed). From 011: H-017 at n = 9 depth 6 negative (E-069).

## What I tried
- 022: stripped the kernel computation from the 019/021 scripts, tallied shape tests on every step; four parallel 480 s runs on 4 cores.
  Wall-clock caps make counts differ between runs (rejecting 1 842 here vs 1 123 in E-097): quote rates, not counts. Prefer `--maxexp`.
- 019: copied scholar's script, added kerdims and a parent tally. 017: monkeypatch wrapper. 014: `batch.py orbits 17`. 011: shards by `xargs -P 4`.

## What I would do next
1. Look for a long-square tilting step off the walks (E-078 `long`, 3-path relations, non-guard steps): that decides whether the shape is a criterion or a walk artifact.
2. n = 6 classes 2-3, n = 7 classes 1+ with `--maxexp`; shape of the 2 n = 8 c2 rejecting parents.
3. Overnight: n = 8 c2 depth 8 to completion (`--budget-sec 1500`), then depth 9; H-017 depth 7 at n = 9.
4. Row-set comparison of n = 15 big orbit; 4-letter scan n = 16 by shards; n = 17 key-coarser lists.
- Watch: a timeout/budget is not a verdict; parallel timings are inflated; equal sizes are not equal sets; hand-built family (E-078) is not the walk's distribution; a loose shape test passing everywhere is weak evidence.
