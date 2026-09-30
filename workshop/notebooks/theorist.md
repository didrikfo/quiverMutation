# Theorist notebook (rewritten round 009)

## What I believe now
- Drift family `aax` (a=3..6), chain c = x+o conserved (E-065 for a=3; 44x/55x/66x same local neighbour pattern at n=20). Rigid (chain pair {c, n-c}, k = 2x+3-a)
  for a = 3, 5, 6 (55x, 66x predicted before running, n = 14, all orbits closed); a = 4 merges: `444 -> 34` (double mutation collapse) and `34` reaches the
  slider `44` and the big orbit O* (3767 at n=14). Criterion = "interior collapse product of the seed is 34" (computed, not derived).
- End link is uniform for a=3,5,6: O* holds exactly `aaa@0` and `aax@hi`; every other orbit is a pure chain pair.
- Weakest step: the path 34 -> 403 -> 34@+2 -> 3333 -> 3403 -> 44 (free/reduced moves). Why only 34? unexplained.
- "Orbit holds some word at all offsets" does NOT predict merging (10/24): orbits of 334@1 hold 36, 66 everywhere and stay rigid. Do not retry.
- Older: k(33x)=2x explained in round 006 (double mutation, not rule table); general k = 2x + w0 - x0; 34x k = x+3 (only 345 explained); 45x no drift;
  H-021' stands; "s in cons" is tautological.

## What I tried
- Round 009: theorist_{nbrs,translate,closure,translator,path34}.py in rounds/009. Closure test n=14 only (55x took ~4 min).
- Not tried: n=15,16 for 55x/66x, 77x, x >= 9, proving the 34 link, 34x by the same collapse lens (344 -> 24, 355 -> 25, 366 -> 26 all non-34).

## Next
- Hand-derive why `34` slides to `44` (double mutation on 3333/3403), and why 23/45/56 do not: that would make the criterion a lemma.
- 55x at n=15, 16 (plan first); 77x at n=15.
- Use the collapse lens for 34x/35x/36x: is `344 -> 24` the reason 344 is tiny (47-118 rows)?
- Blind spots: four values of a with one exception is thin; the criterion is fitted to a=4 by construction until the 34 step is explained.
