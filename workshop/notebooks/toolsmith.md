# Toolsmith notebook

## Rounds 003-102 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s (use `timeout 115 sh -c 'while pgrep ...'` loops), `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); Hom+Cartan has no NO power for iso, the iso meet is the check.
- T10 n = 9 (E-234..E-268): key-kept J != 0 steps exist; join searches cannot test the J = 0 premise. Tools: `087/toolsmith_d6scan.py`, `098/toolsmith_j98.py`, `102/toolsmith_litex.py`.

## Rounds 105-111
- "A_n" = Nakayama rad^(n-1)=0, "B_n" = lone 3-arrow relation; both derived equivalent to kD_n. `107/toolsmith_cones.py`, `107/toolsmith_tiltsearch.py`, `050/skeptic_tilt.py` (independent Hom test). R109: on every shortest B_n -> A_n path gate admits, key kept.
- R111: all 89 D = 0 steps Z-congruent, but Cartan test has no NO power at n = 7.

## Rounds 113-114 (n = 7 class 1, E-201's 13 steps)
- Believe: all 13 steps have Hom-verified tilting chains (6-12 edges; K0 det +-1, loop-free, labelled Cartan = dim Hom, meet = structure-graph iso up to coefficients, generation assumed). Step 8 (R114): `R4 F1 R7 F1 F4 F3 | F1 F3 R7 R4 R3 F5`, needed 6+6 levels (balls ~17k each).
- Why 113 stalled: pure throughput. Fix = `114/toolsmith_s8par.py` (4-core pool, expand smaller ball, resumable /tmp/w/s8_I.pkl); reusable for any step or other class.
- Control power: meet iso does nearly all rejecting (Hom 1/77, K0/Cartan/loop 0); end swaps 12/12 rejected. No NO case at n <= 7, so this is only YES-side consistency; R-candidate not moved.
- Traps: /tmp checkpoints vanish between sessions (069 ckpt needed 2 x 9 min to rebuild); verify on parent-side algebras with triple arrows is slow (A-side corruption control > 10 min).

## Next
- A-side corruption control for step 8 (split per edge into parallel jobs); class 2 record 6 (E-212) via s8par + K0.
- A NO-capable invariant is the real blocker for T10: find any key-equal pair of inequivalent algebras (n = 8, 9?) and test the walk there; key-refusal control on E-289 pipeline; n = 8 cones sizing.
