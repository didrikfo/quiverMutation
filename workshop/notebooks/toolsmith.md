# Toolsmith notebook

## Rounds 003-090 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s, `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings (two parallel builds ran ~1.7x slower), check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); meet up to sign-aware iso. Hom+Cartan has no NO power, the iso meet is the check.
- T10 n = 9: key-kept perI J != 0 steps exist in both F-010 balls and rejoin their own class by J = 0 moves (E-234, E-239, E-244). join97 README misread: None-key drops 9 / 11 (E-250). The 7 + 7 NO says nothing on the J = 0 premise.
- Tools: `087/toolsmith_d6scan.py` (build / parents / scan, `--budget-hours 0.14`, env DEPTH), `090/toolsmith_nonekey.py`, `093/toolsmith_unionjoin.py` (lna | extra | meet | --plan), `098/toolsmith_j98.py` (jedges | meet | control).

## Rounds 093, 095, 098
- 093 (E-255): 2400230(8) vs 3304000(6), F/R + lnaMoves layer: MISS; 2400230 isolated in lnaMoves. 095: re-issued, I only planned (lesson: read notebook + STATE T10 line first).
- 098: ran it for real (rebuilds ~13 min wall, 2 shards each). All 8 key-kept J != 0 edges (1 from A depth 8, 7 from B depth 5-6) as undirected moves + radius-1 balls: 0 WL-bucket pairs, no meet; all 8 children non-LNA. Control: 4 of 7 B children re-join B at 6, rest at 7-8 (E-244). Verdict: move-set statement; nothing on the J = 0 premise.
- Believe: the 2400230 vs 3304000 calibration is exhausted at feasible radii; larger radius cannot help (balls share no WL hash). Not an OVERNIGHT case.
- Traps: scripts exec the 087 prelude (guard `A_`, `P`); `background &` subshell makes the task notice fire at once, poll the logs; A8 build needs a rerun after exit 2.

## Next
- Only worth doing: a non-move-set certificate for a non-LNA pair (theorist/skeptic), or `lnaMoves.discoverMoves` with a seed row (2400230) explaining why no verified move exists (small, with test).
- 3304000 depth-7 / 2400230 depth 9-10 scans only widen a MISS; do them only if asked.
