# Maverick's notebook

## What I now believe (after round 027)
- S-1 first sitting (n = 8, 9, 10): deletion keeps same-class pairs together about 2x chance (0.43 vs 0.19 at n = 9, i = j = middle) and no position (i, j) works; a "mid + fewest covering relations" rule lands 78% of LNAs in their class's modal image.
- What does work: delete a free vertex at an end (head or tail run) when the run has >= 3 free vertices: image class is a function of source class at n = 8, 9, 10 (3/3, 5/5, 10/10 classes). K = 2 fails once or twice, K = 1 often. Not derived; may be implied by H-020/free move.
- The "some deletion agrees" coverage bound is trivially ~1 (99.8%): do not quote it.
- Class labels: key, cospectral keys split by orbit+mirror; 16 LNAs at n = 9 (176 at n = 10) stay unresolved in the cospectral keys and are dropped.
- Earlier (r023) H-017 beliefs stand: relations > cords among quipus survive depth 6 at n = 9; Euler signature pos(C+C^T) <= n-2 iff outside every quipu class; mirror-chain depth rule is a fit (0 mismatches n = 8..10, 8/8 out of sample).

## What I tried
- r027 scripts `rounds/027/maverick_{classes,delete,null,rules,dist,free,freetype,endstrip}.py`; all under 1 min.
- r018-r023: cord producer, `maverick_predict*.py`, `maverick_chain.py`, `maverick_twobig.py`.

## Watch for
- Pair statistics are dominated by the big classes; compare with the all-pairs null.
- `freeMoves.derivedOrbits(rules=None)` still leaves cospectral orbits unmerged: do not trust orbit identity as class at n >= 9 in cospectral keys.
- A negative at depth d excludes only distance <= d.

## Next
- Resolve cospectral unresolved LNAs, run n = 11 for the free-end K >= 3 statement; derive why K >= 3 (and why K = 2 fails for the 300-member n = 9 class).
- Question 3 of S-1: LNAs whose core sits with room (head >= 3) vs the same core at head < 3: is class equal across lengths when head >= 3 (a "stable core" statement)?
- n = 12 two-big LNAs for the mirror chain (carried); depth 7 at n = 9 overnight for H-017.
