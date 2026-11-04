# Experimentalist notebook (rewritten round 114)

## What I now believe
- S-1/H-026: K-threshold K0 = (n-1)/2 - 2 holds at key level for the lone-3 class at n = 17, 19, 21, head (113) AND tail deleted directly (114, predicted first, 0 mismatches vs pwh.removeVertex(n,r,n)); tables identical up to key labels. Scope: not class level, not tested below K0-1 at n = 19, 21.
- Reply 114: E-263 (4,7) "12 keys" = writer bug (''.join of entries >= 10 -> 440 malformed lines in h=4 hit files; (5,6) h=5 files clean). (4,7) has K0 = 4. Non-central lone-3 classes: K0 fits min(h, n-4-h) (n=13 h=2,3,4; n=15 h=4,5; n=17 h=5,6). H-026 wording falsified; (n-1)/2-2 only for the central class. Fit on 7 points, untested at n = 19 h<7.
- Counts 780/265/129/65 (+ 2^m style) recur across classes; unexplained.
- T10/H-015 (110): key guard refuses real steps at n = 6 class 0 (1150 J != 0 steps), none at J = 0; no guarded walk known to leave a class. Open: n=7 c1, c2 record 6.
- T7/H-010 (107): no lowering control for start-gap-2 pair at k <= 4. T6 (105): defect not constant per key class at n = 9 depth 5.
- Practical: tail row = shifted row + zeros at the end, image = drop last entry. Write outputs to files; log predictions in a file first.

## What I tried
- 114: `rounds/114/experimentalist_tail.py` (TAIL=1, h=). 113: `_kthresh.py`. 110: `_refusal.py`. 107: `_chain.py`; 105: `_defect_census.py`.

## What I would do next
1. Test K0 = min(h, n-4-h) at n = 19 (h=5,6), n = 15 h=3; n = 23 tail.
2. T10: n = 7 class 1 Hom-verified chain (overnight). T7 chain at k=5; T6 depth-6.
