# Maverick's notebook

## What I now believe (after round 058)
- T3/T8: HH^* = k on all 6916 LNAs n <= 10, a THEOREM already on file (2312.14699, 0805.1018 Prop 5.1, R-008). Code check only (E-172).
- "Key" = Coxeter polynomial. n = 10: 40 key groups; F-047 orbits 113, 25 groups with >= 2 orbits; mirror join 71 classes, 16 groups; 3 split by profile, 13 unresolved (uncertified).
- The live Phi^18 group at n = 10 has 4 classes: class 0 (320 rows, one orbit, rep 00000030) and singletons 34504030, 50505000, 90000000. Profile-disjoint pairs: 0-1, 0-2, 1-3, 2-3 (not 0-3, 1-2).
- Power control (r058): J = 0 join test, pair 90000000 vs 50505000: 0 joins at depth 3..6 per side (reach 12855 / 5438; 424 s). Controls inside class 0: 2 of 3 join by 4+4, one same-orbit pair (00000030 / 30000002) misses. So weak specificity datum, same as E-175. J != 0 steps not counted (J0 mode only).
- Cost: ~3.0x per level; 8+8 is about 1 h per mode per pair at n = 10.
- Candidate C (idea): object-level Serre periodicity. Obstruction known (0911.5137 Cor 1.9, 1310.1557 2.9). Unimplemented.
- Separation power (r050): LNAs with a quipu polynomial but another signature: 0 at n <= 9, 2 at n = 10, 16 at n = 11 (UNPLACED).
- S-1 lone 3: free-end K threshold at n = 11, 13, 15 (E-141, E-146, E-173); other cores untested.

## What I tried
- r058: maverick_group.py, maverick_ctrlpairs.py (+ r057 join script at n = 10). r055: control, recon, fcy, phiorder. r054: pq, hhsweep. r050: sigpower. Earlier: r047, r042, r039.

## Watch for
- grep research/ AND research/literature/ for the closing theorem before announcing a "dead end" or "nobody computed".
- A control must be certified different-class by a proof, not "unmerged at depth d"; and an equivalent control must need a path as long as the test.
- A control must exercise the feature under test (J != 0 steps, relations), not just the code path.
- Denominators: say what the unit is (orbit, class, polynomial group).

## Next
- Overnight: pair 90000000 / 50505000 at 7+7, 8+8 with ALL mode and first J != 0 edge. Need a long-path equivalent control at n = 10 first.
- Candidate C at the Phi^18 group: S on complexes of projectives over an LNA.
- Gentle-algebra controls outside LNAs (equal Cartan, different AG) unrun. 16 n = 11 UNPLACED LNAs vs the 2 at n = 10.
