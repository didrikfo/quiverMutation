# Maverick's notebook

## What I now believe (after round 042)
- S-1 lone 3, key level: the free-end K-threshold law now has three lengths. n = 11: K >= 4 holds, K >= 3 fails (E-118/E-133). n = 12: K >= 4 holds, K >= 3 fails (r042; class of (3,5)/(5,3), 2746 LNAs, ONE orbit; K >= 3 ends go to 2 keys, K >= 4 ends to 1). n = 13: K >= 4 fails (E-139; two orbits 4349 + 674). Single-3 key table: K0 = 3 first fails at n = 11, K0 = 4 at 13, K0 = 5 at 15 (unrun at class level).
- So the K >= 3 failure is not an orbit-structure effect (one orbit at n = 12); it is the lone 3 itself with h != K. The self-mirror lone 3 at (4,4) never fails.
- The failing witness is always the lone 3 (two ends, different keys); the class-level content is that nothing else in the class breaks K >= 4 at n = 12.
- Compatible deletion rule for lone relations: delete from the shorter free side (mirror-aware); untested for other cores.
- Earlier: n = 11 failing class (1305 LNAs) I1/I2 split = lone 3 or 7 with K_eff = 3; mirrorRow redo 94 ends, 47 keys, 0 conflicts; E-125's 82 vs my 94 unexplained. H-017 beliefs of r023 stand.

## What I tried
- r042: maverick_n12.py (scan 2 shards 1 min, ends+orbits 1 min); fixed maverick_single.py label block (E-115 label dicts are keyed by length n-2 rows). Full n = 12 scan is cheap (58786 rows).
- r039: n = 13 class, ends, orbits. r037: mirror redo, normal form, single-relation keys. r034, r030, r027, r018-r023 earlier.

## Watch for
- Image comparisons are by key, not label (no E-115 labels at n = 11 for these classes beyond the table). "Differs" sound, "same" not.
- Orbit within a key class = lower bound for the derived class.
- One class per length is a small sample; the lone 3 is the simplest core, not typical.

## Next
- Non-lone cores (pairs of relations, lone 4 / 7): does K0 shift? K = 2 failures at n = 11 (9 classes).
- n = 15 for K0 = 5: build the class from the lone-3 orbit, not a scan (n = 13 scan was 4 x 5 min); size first.
- Depth 7 at n = 9 for H-017; n = 12 two-big LNAs (carried).
