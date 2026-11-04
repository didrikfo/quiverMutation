# Toolsmith notebook

## Rounds 003-102 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s, `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); Hom+Cartan has no NO power for iso, the iso meet is the check.
- T10 n = 9 (E-234..E-268): key-kept J != 0 steps exist; join searches cannot test the J = 0 premise. Tools: `087/toolsmith_d6scan.py`, `098/toolsmith_j98.py`, `102/toolsmith_litex.py`.

## Rounds 105-109
- "A_n" = Nakayama rad^(n-1)=0, "B_n" = lone 3-arrow relation; both derived equivalent to kD_n. `107/toolsmith_cones.py` cone closure (n=7: 2m38s); `107/toolsmith_tiltsearch.py` BFS over tilting exchanges. `050/skeptic_tilt.py` = independent chain-map Hom test.
- R109: on every shortest B_n -> A_n tilting path the gate admits and the key is kept (`109/toolsmith_pathaudit.py`); match there is by labelled Cartan only.

## Round 111 (n = 7 class 1, E-201's 13 steps, Cartan route)
- All 89 D = 0 steps Z-congruent, but the test has no NO power at n = 7 (`111/toolsmith_c1_13_exact.py`, pari-free). Box search over small combinations misses real congruences; a miss is not NO.

## Round 113 (Hom-verified chains for the 13 steps)
- Believe: 12 of 13 steps (all but 8) have explicit Hom-verified tilting chains (6-11 edges, K0 det +-1 on every edge, labelled Cartan = dim Hom, meet = structure-graph iso). Steps 3, 9 needed depth 5+6 (LIM=5,6, cap 150000). Step 8 (parallel arrows): slow, no meet at 5.5 k/9.8 k algebras after ~4 resumes; undecided.
- Tool: `113/toolsmith_c1_chain.py list | chain I depth cap | ctl I` (resumable ckpt /tmp/c1chain_ckI.pkl, env LIM=a,b per-side depth); `113/toolsmith_c1_k0ctl.py`. K0 check added (`k0det`); duplicated summand det 0 in 708/708, real +-1 in 118/118.
- Controls: single-edge corruptions rejected except 2/step that are genuine alternative edges; end swapped for another step's child rejected 36/36. No power for "leaves the class".
- Traps: checkpoints and logs in /tmp survive container restarts only sometimes (the earlier run's logs were in /tmp/claude-0); tool calls cut at 120 s even with `timeout 590`, so poll in <= 115 s loops; `ls step_*` sorts 10 before 2.
- Meet iso ignores relation coefficients (arrow identities collapsed in relation terms).

## Next
- Step 8: deeper search (LIM=6,7) or cone closure from the child side; class 2 record 6 (E-212) through the same verifier with K0.
- Key-refusal control on E-289 pipeline; n = 8 cones sizing.
