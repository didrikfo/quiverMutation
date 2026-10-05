# Maverick's notebook

## What I now believe (after round 034)
- S-1 free-end deletion: with E-115 labels, "K >= 3 free vertices at the end => image class is a function of source class" holds n = 8, 9, 10 (3/3, 6/6, 12/12) and fails at n = 11 in ONE class (key (1,1,0,-1,-2,-3,-3,-2,-1,0,1,1), 1305 LNAs, orbits 15107 (943), 15035 (362)). K >= 4 holds (10/10 resolved at n = 11). K = 2 fails 0/1/2/9 classes at n = 8..11.
- r034 per-end table (`rounds/034/maverick_endtable.py`, 80 s): 82 ends with K >= 3 (66 at K = 3, 16 at K = 4), head/tail mirror-identical. Image is a function of (K, oriented core word): 77 keys, 0 conflicts. 12 core words occur at K = 3 and K = 4: K = 3 -> image I1 (n = 10 key ...-2,-2,-2...), K = 4 -> I2 (...-2,-3,-2...). So the failure is "room to move": deleting from a run of 3 leaves a run of 2 (the K = 2 regime), from 4 leaves 3. Not a label artefact (image keys differ; both images inside orbit 15107). Words ending in 3 go to I1 at K = 3 (25 words), others to I2.
- Earlier: pair statistics of deletion 2x chance, no positional rule (E-112); H-017 beliefs of r023 stand (relations > cords among quipus to depth 6 at n = 9; Euler signature pos <= n-2 iff outside every quipu class; mirror-chain depth rule is a fit).

## What I tried
- r034 endtable (above). r030 `maverick_endstrip2.py N`, `maverick_fail11.py`. r027 `maverick_*.py`; r018-r023 cord producer, predict, chain, twobig.

## Watch for
- Pair statistics are dominated by big classes; compare with all-pairs null.
- Key-only class labels can merge derived classes (source side); orbit-split images inside one orbit is the clean test.
- Head-end orientation in r034 was a digit-string reversal, not freeMoves.mirrorRow; rerun if it matters.
- A negative at depth d excludes only distance <= d. One failing class is a small sample for a "law".

## Next
- Check the same K-vs-(K-1) same-core pattern for the 9 K = 2 failing classes at n = 11.
- n = 12 (orbit-only test): does the K = 3 failure persist for the same core words, or does the threshold rise? Needs n = 12 free-move orbits; size with --plan.
- Derive why last letter 3 matters (H-020 rule table; theorist).
- S-1 Q3 in the other direction (simple LNA reached only via room) untouched.
- n = 12 two-big LNAs for the mirror chain (carried); depth 7 at n = 9 overnight for H-017.
