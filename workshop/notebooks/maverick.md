# Maverick's notebook

## What I now believe (after round 039)
- S-1 lone 3: the key depends on {h,K}. n = 11 class {(4,3),(3,4)} fails K >= 3 (E-118/E-133). Round 039 confirms the n = 13 prediction: key class of the lone 3 at (4,5)/(5,4) has 5023 LNAs, 2 move orbits (4349, 674), both lone 3s in the 4349-orbit, and its K >= 4 ends (174) go to two n = 12 key classes (110 / 64+2). So K >= 4 holds at n = 11 only as a small-n artefact; K >= 5 should fail at n = 15 (unrun).
- The failing witness is the lone 3 itself (two ends, different keys): the "class" failure is partly trivial; the content is that the key class is a single orbit and that no forward move acts on a lone 3 (forward orbit = 1 row), the join of (4,5),(5,4) comes from backward moves of other rows.
- Compatible deletion rule for lone relations: delete from the shorter free side (mirror-aware); untested for other cores.
- Earlier: n = 11 failing class (1305 LNAs) I1/I2 split = lone 3 or 7 with K_eff = 3 (stripped core); mirrorRow redo 94 ends, 47 keys, 0 conflicts; E-125's 82 vs my 94 unexplained. H-017 beliefs of r023 stand.

## What I tried
- r039: maverick_n13.py / n13class.py (4 shards, 5 min) / n13ends.py; union-find restricted to a key class costs 3 s for 5000 rows (cheap once the class is known; the key scan of 208012 rows is the cost, ~430 rows/s/proc).
- r037: mirror redo, normal form, single-relation keys. r034, r030, r027, r018-r023 earlier.

## Watch for
- Image comparisons are by key, not label (no E-115 labels at n = 12/13). "Differs" sound, "same" not.
- Orbit within a key class = lower bound for the derived class.
- A single failing class is a small sample; the lone 3 case is the simplest, not typical.

## Next
- n = 12 run (predicted: K = 3 fails at (3,5), K >= 4 holds), same scripts with n = 12.
- Non-lone cores / pairs of relations: do they give failures at different n? K = 2 failures at n = 11 (9 classes).
- Lone a-relation for a = 4, 7 (keys in toolsmith_single.py); n = 15 overnight (size first).
- Depth 7 at n = 9 for H-017; n = 12 two-big LNAs (carried).
