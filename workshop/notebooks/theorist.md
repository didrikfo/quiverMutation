# Theorist notebook (rewritten round 017)

## What I believe now
- Round 017: orbits of this move set are labelled by a class J = offsets of `333@o'` held (pairs {j, n-6-j}; E-065's drift label
  j = x+o-3 for `33x@o`). `3x` has j = x-4 (x = 4..7, all placements at odd n; fold at x = 8, 9). `444` and `34` are j = 0 (the
  "444 orbit"), `35`, `55`, `455`, `2455`, `3334` are j = 1 (the small `235/255/455` orbit). Checked n = 12..16.
- Lemma R (double mutation at the middle of `a b b d`): `(a,b,b,d) -> (a-1,b,d+1)`, `(a,b,b) -> (a-1,b)` for 3 <= a <= b: tested 250/250 and
  119/119 at n = 14. So `444 -> 34` (j 0) but `3334 -> 35` (j 1): the 4 is covered by `r` and becomes a 5. "Letter 4" vs "collapse to 34"
  is settled for these words: it is neither; it is the tail after R.
- `k(33x) = 2x` follows from the class table (`o + o' = n - 2x`); E-065's upper bound is enumerated n = 12..16, still not proved.
  No invariant known that separates classes (nulls from 013 stand): "class" = closed orbit under this move set, not derived inequivalence.
- Round 015 (T5): the n = 8 class 2 "gate True, tiltingPlus True" steps are a defect of `arrowPaths.reduceAgainstPivots` (not a normal
  form; step 7 drops a relation); monkeypatch in `rounds/015/theorist_fix.py`; toolsmith owns the library fix. Cartan congruence agrees with
  `tiltingPlus` on 11 replayed rejections.
- Earlier: P (3@2/3@3) and Q (5@0/6@0) each closed under shift by 2 (013); P and Q are the parity pair at even n and merge at odd n
  (017 tables: `35` at odd n one orbit, at even n two).

## What I tried
- 017: `theorist_{small,single,chain,label,class,rrule}.py`. Fastest: label every placement by (orbit size, J) and print the table; labelled
  `theorist_path.bfs` gives the move names; hand-apply the double mutation to intervals to read off the lemma.
- 015: replay of recorded walks, wrapped `_kernelOverIdeal`. 013: labelled BFS.
- Nulls (do not retry): no GF(2) functional, integer statistic or SNF separates P from Q.

## Next
- Words with J: predict J for every 4-letter word with a 4 by (apply R, read `3x -> x-4`) and compare with the scans at n = 13, 15.
- Write `3x@0 -> 33(x-1)@0` (anchored) as a rule row; derive J for `3x` for all x. Explain the fold and the x >= 8 range.
- Why is the shadow orbit (no `333`) a parity class (even n: 272 / 446 / 148 sizes)? Same question as P/Q.
- Library: after the `reduceAgainstPivots` fix, ask whether `tiltingPlus`'s rank step is affected.
- Blind spots: J is a name for an orbit under a finite move set; I have no inequivalence proof. Lemma R tested at n = 12, 14 only for
  interior placements; boundary-touching cases are not covered.
