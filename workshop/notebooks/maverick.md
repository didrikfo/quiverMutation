# Maverick's notebook

## What I now believe (after round 121)
- n = 12 key group of 47 (key (1,1,-2,-4,-2,3,6,3,-2,-4,-2,1,1)): SNF(S) (0,1^10,20) on 14, (0,1^10,12) on 33, now double-checked (flint + pari). The move census (`derivedOrbits(12, free, edges, doubles)`) gives 9 orbits: 13, 13, 7 (33-class), 4, 4, 2, 2, 1, 1 (14-class); 0 mixed, 0 moves leaving the group; 136 in-group census edges all keep SNF. So: no census merge across the split, no bug.
- Key-guarded balls (library DFS, both orientations): 3 roots in 14-class to depth 4-5 (one to depth 6 right side), 2 roots in 33-class to depth 4: every node keeps root SNF. Zero crossings.
- BUT the guard never fires here: guard-off ball = guard-on ball (depth 3-4, two roots). So these balls test the gate, not the guard, and there is no walk-level NO control (no step in reach is known to leave the class). Detector power is planted only (33/33 flagged, 0/14 false).
- Balls reach almost no other LNA (root, opposite, free variants). LNA-level "lands among the 33" is a weak target; node-level SNF is the real detector.
- Within each SNF class the census leaves 6 and 3 orbits: H-032's second falsifier (walk components the Cartan level does not separate) lives there.
- E-269 / H-024 (r119): NO power at one n = 11 group (lam 1.458, det S 8); no planted pairs at n = 12 low-lam groups by entry perturbation; singular search YES-only (no bound).
- S-1 / H-026 (r115): K0 = min(h, n-4-h) at odd n; even n fails at the self-mirror class (E-307).
- H-017b text written (r119), overnight proposal, not run.

## What I tried
- r121: workshop/rounds/121/maverick_census47.py (84 s; CHECK=1 whole-length 111 s), maverick_ball47.py (ROOT, DEPTH, SIDES, CTRL=off|planted; per side 45 s at d4, 145 s at d5, ~440 s at d6, plus 84 s lnaKeys table).
- r119 plant2/plant3/qsgroups/sing/sing_cross; r115 mistat; r114 defectstat; r111 k17; r109 planted; r107 n12low; r105 lowlam; r099 iso140; r097 cartanZ.

## Watch for
- Venv loses cypari/python-flint: `uv pip install --python .venv/bin/python cypari python-flint` (package is `cypari`, not `cypari2`). flint fmpz_mat has .snf().
- lnaKeys(12) costs 84 s per process; cache members if running many balls.
- 10 min per command; background + until-loop to wait (plain sleep blocked).
- A NO control must have an independent certificate; a guard-off ball that never moves the key is NOT a control.
- Opposite-side LNAs in ball47 are reported as found on the dual (mirror); labels are mirror-invariant in this group.

## Next
- Find a root (any n) where the unguarded walk DOES move the key, to get a walk-level NO control; or graft a planted crossing step (toolsmith).
- H-032 inside one SNF class: meetingPoints / E-314 singular search between the 6 orbits of the 14-class and the 3 of the 33-class.
- Overnight if approved: depth 6-7 both sides from 3040500030, 0040033030.
- Literature slot (deferred from r121): Ladkani, BGS/Avella-Alaminos-Geiss on complete invariants.
- Even n S-1; H-017b census when placed.
