# Experimentalist notebook (rewritten round 114)

## What I now believe
- S-1/H-026: K-threshold K0 = (n-1)/2 - 2 holds at key level for the lone-3 class at n = 17, 19, 21, head (113) AND tail deleted directly (114, predicted first, 0 mismatches vs pwh.removeVertex(n,r,n)); tables identical up to key labels. Scope: not class level, not tested below K0-1 at n = 19, 21.
- n = 15 class (4,7)/(7,4) (114): 3800 rows at K >= 2, head = tail, splits for K = 2..4, one key K >= 5 (K0 = 4). E-263 says this class is "messy, 12 keys at K >= 5": unreconciled (different K bookkeeping / sub-family?). Do not claim a threshold for it until reconciled.
- Counts 780/265/129/65 (+ 2^m style) recur across classes; unexplained.
- T10/H-015 (110): key guard refuses real steps at n = 6 class 0 (1150 J != 0 steps), none at J = 0; no guarded walk known to leave a class. Open: n=7 c1, c2 record 6.
- T7/H-010 (107): no lowering control for start-gap-2 pair at k <= 4. T6 (105): defect not constant per key class at n = 9 depth 5.
- Practical: tail row = shifted row + zeros at the end, image = drop last entry. Write outputs to files; log predictions in a file first.

## What I tried
- 114: `rounds/114/experimentalist_tail.py` (TAIL=1, h=). 113: `_kthresh.py`. 110: `_refusal.py`. 107: `_chain.py`; 105: `_defect_census.py`.

## What I would do next
1. Reconcile (4,7) at n = 15 with E-263 (run K1=5 over all rows, compare K definitions); n = 23 tail.
2. T10: n = 7 class 1 Hom-verified chain (overnight). T7 chain at k=5; T6 depth-6.
