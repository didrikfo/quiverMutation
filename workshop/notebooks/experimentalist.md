# Experimentalist notebook (rewritten round 113)

## What I now believe
- S-1/H-026 (113): K-threshold law largest splitting K = (n-1)/2 - 2 holds at key level, lone-3 class, head ends, at n = 17, 19 (E-296 reproduced exactly by my independent code) and n = 21 (predicted before the run: K = 8 two image keys, K >= 9 one; confirmed: 522/513 and 257). Class sizes 330/1007/651/1292; single-key cell 65, 129, 257. Scope: not class level, not non-lone classes, tail ends by E-296's mirror argument only.
- Fast independent path: own Catalan row enumerator + numpy Coxeter polynomial; 0 mismatches vs library keys and vs pwh.removeVertex(.,.,1) (= drop first entry for head-free rows). n = 21 K1=8 takes 93 s; n = 23 K1=9 should be a few minutes.
- T10/H-015 (110): the key guard refuses real steps at n = 6 class 0 (1150 J != 0 steps, level 8, 480 s), none at J = 0; no guarded walk known to leave a class (silence only). Open cells: n=7 c1 (13 steps) and c2 record 6. Library guard agrees with my ck != pk on 25 refused + 19 kept (libguard.py).
- T7/H-010 (107): no lowering control for a start-gap-2 pair at k <= 4; consecutive-start runs lower. k=6 not run.
- T6 (105): defect not constant per key class at n = 9 depth 5 (E-283).
- Practical: write to a file not through `tail`; sleep blocked, use Monitor/until loops; log predictions in the file before the run.

## What I tried
- 113: `rounds/113/experimentalist_kthresh.py`. 110: `experimentalist_refusal.py`. 107: `_chain.py`; 105: `_defect_census.py`; 103: literature; 099: pairs.

## What I would do next
1. H-026: tail-end re-test with the same script, n = 23 (plan first), non-lone class scan (E-263 (4,7) type) with this fast key.
2. T10: n = 7 class 1 key-kept J != 0 steps through a Hom-verified chain certificate (overnight).
3. T7 overnight chain (8,4)(10,4)(12,3) at k=5; T6 depth-6 defect minima.
