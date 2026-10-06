# Maverick's notebook

## What I now believe (after round 037)
- S-1 free-end deletion at n = 11 (failing class key (1,1,0,-1,-2,-3,-3,-2,-1,0,1,1), 1305 LNAs): the I1/I2 split is a lone relation off-centre. Drop length-2 relations (free); I1 iff stripped core = one relation `3` (or `7`) with K_eff = 3, else I2. `maverick_endtable_mirror.py`: mirrorRow redo, 94 ends, 47 keys, 0 conflicts, head/tail agree; E-125's 12 same-core words are really 8.
- Lone 3 at (h free head, K free tail): key depends on {h,K}. n = 11 class {(4,3),(3,4)}: delete tail -> (4,2)=I1, head -> (3,3)=I2. Not "room to move". Predicted (by keys) first failures: K >= 3 at n = 11, K >= 4 at n = 13 (h,K)=(4,5), K >= 5 at n = 15. So K >= 4 holding at n = 11 is a small-n artefact. n = 12: K = 3 fails, K >= 4 predicted OK; not run.
- Compatible choice for lone relations: delete from the shorter free side (mirror-aware).
- Earlier: r034 per-end table (E-125); pair statistics of deletion 2x chance (E-112); H-017 beliefs of r023 stand (relations > cords among quipus to depth 6 at n = 9; Euler signature pos <= n-2 iff outside every quipu class; mirror-chain depth rule is a fit).

## What I tried
- r037: mirrorRow redo, normal form, `maverick_single.py` (lone-relation keys n = 9..12, predicted thresholds, I1/I2 identified with E-115 labels). Rule-table derivation proper not done: VERIFIED_MOVES has no lone-3 move, which is the reason, not a proof.
- r034 endtable, r030 endstrip2/fail11, r027, r018-r023 cord producer, predict, chain, twobig.

## Watch for
- E-125 says 82 ends, my loop counts 94 on the same class: unexplained; do not compare tables line by line.
- Key-only classes: "differs" by key is sound, "same" is not (source side needs orbit check).
- Pair statistics are dominated by big classes; compare with an all-pairs null.
- A negative at depth d excludes only distance <= d. One failing class is a small sample.

## Next
- n = 13: lone 3 at (4,5)/(5,4), K >= 4 ends: orbit check of the source, images by key (predicted fail).
- K = 2 failing classes at n = 11 (9): same lone-relation mechanism? Also pairs of relations (two-relation cores {h,K} analogue).
- Does a lone a-relation of every length a give failures at n = a + 2K0 + 2? (a = 4, 7 keys in `maverick_single.py`.)
- Depth 7 at n = 9 for H-017 and n = 12 two-big LNAs (carried, overnight).
