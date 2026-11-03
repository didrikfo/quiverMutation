# Experimentalist notebook (rewritten round 110)

## What I now believe
- T10/H-015 (110): the key guard refuses real steps: n = 6 class 0 (2 LNAs) refuses 1150 J != 0 steps (first at BFS level 8, 480 s), refuses no J = 0 step; classes 1-3 and n = 5 show no refusal and no key-kept J != 0 step. Key-kept J != 0 steps first occur at n = 7 c1, c2 (E-147); 8 of 9 c2 ones are in-class (E-212). So no guarded walk is known to leave a class (silence only: 5 of 6 cells have no J != 0 step; 'cannot leave its class' was withdrawn after review); the open cells are n=7 c1 (13 steps) and c2 record 6.
- Library guard (`mutationSearchDepthFirst` depth 1, guard on/off) agrees with my `ck != pk` on 25 refused + 19 kept sampled steps (libguard.py). 480 s counts are load-dependent.
- Refusal needs a long square; J != 0 is rare under level 8, so shallow audits (depth 6) see "guard never fires".
- T7/H-010 (107): no lowering control for a start-gap-2 pair at k <= 4; consecutive-start runs lower; (8,10,12) chains stay at base. k=6 not run.
- T6 (105): defect not constant per key class at n = 9 depth 5 (E-283).
- Practical: write to a file, not through `tail`; the sleep command is blocked, use Monitor or `until` loops in background; key classes at n=6 number 4.

## What I tried
- 110: `rounds/110/experimentalist_refusal.py` (classes 0-3 at n=6, 0-1 at n=5). Reused `001/scholar_h015.py`, `041/experimentalist_keyoff.py` design.
- 107: `experimentalist_chain.py`; 105: `_defect_census.py`; 103: literature scripts; 099: pairs.

## What I would do next
1. T10: n = 7 class 1 key-kept J != 0 steps through a Hom-verified chain certificate (overnight; E-201 walk 27 min).
2. T7 overnight chain (8,4)(10,4)(12,3) at k=5.
3. T6 depth-6 defect minima; gate check on length-3 relations (Ladkani 1001.4765).
