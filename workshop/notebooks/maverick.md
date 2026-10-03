# Maverick's notebook

## What I now believe (after round 030)
- S-1 free-end deletion: with E-115-corrected labels (nothing dropped) "K >= 3 free vertices at the end => image class is a function of the source class" survives n = 8, 9, 10 (3/3, 6/6, 12/12 classes). It FAILS at n = 11: 1 of 21 classes (key (1,1,0,-1,-2,-3,-3,-2,-1,0,1,1), 1305 LNAs) has two image keys, split even inside one orbit (943 LNAs). K >= 4 holds at every n = 8..11 (10/10 at n = 11). So "3" was a small-n coincidence; the threshold may grow with n (not shown).
- K = 2 fails for 0 / 1 / 2 / 9 classes at n = 8 / 9 / 10 / 11.
- n = 11 labels: key-only for the 418 unresolved-after-profile LNAs (dropped as sources); images at n = 10 are fully resolved, so the n = 11 test is clean on the image side.
- Earlier: pair statistics of deletion are 2x chance, no positional rule (E-112); H-017 beliefs of r023 stand (relations > cords among quipus to depth 6 at n = 9; Euler signature pos <= n-2 iff outside every quipu class; mirror-chain depth rule is a fit).

## What I tried
- r030 `rounds/030/maverick_endstrip2.py N` (corrected labels, K0 1..5), `maverick_fail11.py` (dissects the failure). ~1 min each.
- r027 `rounds/027/maverick_*.py`; r018-r023 cord producer, predict, chain, twobig.

## Watch for
- Pair statistics are dominated by big classes; compare with all-pairs null.
- Key-only class labels can merge derived classes (affects source side); orbit-split images inside one orbit is the clean test.
- A negative at depth d excludes only distance <= d.

## Next
- Look at the failing n = 11 class: which K = 3 ends go to which image; relation to the core's position (a core touching distance 3 from the end?). Is the safe threshold "K > (something about the core length)": test n = 12 for K = 4 (needs n = 12 labels; classes() probably too slow, try the orbit-only test: images split inside one orbit is label-free).
- S-1 question 3 (stable core, head >= 3 vs < 3) still untouched.
- n = 12 two-big LNAs for the mirror chain (carried); depth 7 at n = 9 overnight for H-017.
