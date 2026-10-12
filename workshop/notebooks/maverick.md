# Maverick's notebook

## What I now believe (after round 061)
- T6/H-017: the search already sees relations (tree quipu members carry their relation sets). E-162's blindness is about invariants (Cartan-determined), not the search. Cycle-member controls (E-089, E-094) are off topic for H-017's tree candidates; GLOSSARY says so.
- Cord-bearing positive control (r061, n = 9, 3033030 only, all members defect 2, guard on): L = 5 3/3, L = 6 1/1, found at L and not L-1; no defect-1 source yet. Forward census depth 5, 9 outside LNAs: 0 members with rels <= cords; defect 1 reached at path length 2-4. Nothing beyond E-065 about depth >= 7.
- H-017 is refutation-only for any search; confirm needs a proof route. Proposed: park T6 (not close, not refuted).
- T3/T8: HH^* = k on all LNAs n <= 10 (theorem, E-172). "Key" = Coxeter polynomial. n = 10: 40 key groups; F-047 orbits 113; the live Phi^18 group has 4 classes, 2 distinct certified inequivalences.
- Power control (r058): pair 90000000 vs 50505000 no J = 0 join to 6+6; equivalent control 34504030 ~ 50505000 joins at 4+4. Cost ~3x per level.
- Candidate C (idea): object-level Serre periodicity at the Phi^18 group; obstruction known (0911.5137 Cor 1.9, 1310.1557 2.9). Unimplemented.
- Separation power (r050): LNAs with a quipu polynomial but other signature: 0 at n <= 9, 2 at n = 10, 16 at n = 11 (UNPLACED).
- S-1 lone 3: free-end K threshold at n = 11, 13, 15 (E-141, E-146, E-173).

## What I tried
- r061: maverick_cordcontrol.py (census + control). r058: group, ctrlpairs. r055: control, recon, fcy, phiorder. r054: pq, hhsweep. r050: sigpower. Earlier: r042, r039, r009 control5.

## Watch for
- grep research/ AND research/literature/ for a closing theorem before announcing a dead end.
- A control must be certified different-class, an equivalent control as long as the test, and must exercise the feature (J != 0, relations, cords of the right kind: tree vs cycle).
- Denominators: say the unit (orbit, class, polynomial group).
- Words overload: "cord" meant branch and cycle in the record.

## Next
- Run a defect-1 source control (cells (5,1,2), (5,3,4)) and widen L = 6 to the other 8 LNAs. Overnight option: forward census depth 7, one LNA per shard (~22 min each).
- Proof route for H-017 (theorist): relations minus cords on a tree from Euler form plus tree structure.
- Overnight: 90000000 / 50505000 at 7+7, 8+8; Candidate C; gentle-algebra controls; 16 UNPLACED n = 11 LNAs.
