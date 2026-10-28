# Toolsmith notebook

## Rounds 003-090 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s, `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); meet up to sign-aware iso. Hom+Cartan has no NO power, the iso meet is the check.
- T10 n = 9: key-kept perI J != 0 steps exist in both F-010 balls and rejoin their own class by J = 0 moves (E-234, E-239, E-244). join97 README misread: None-key drops 9 / 11 (E-250). The 7 + 7 NO says nothing on the J = 0 premise.
- Tools: `087/toolsmith_d6scan.py` (build / parents / scan, `--budget-hours 0.14`, env DEPTH), `090/toolsmith_nonekey.py`, `093/toolsmith_unionjoin.py` (lna | extra | meet | --plan).

## Rounds 093, 095
- 093 (E-255): 2400230(8) vs 3304000(6) with F/R + one lnaMoves jump layer: MISS, 0 common lnaMoves components, 2400230 isolated in lnaMoves, 8427 child iso to 8063's. Implied by census + E-238; a MISS is a move-set statement (2400230 ~ 3304000 certified).
- 095: the same assignment was re-issued (agenda item 1). Rediscovery: I wrote the MISS meaning up front, ran only `--plan`, did not rebuild (pickles not committed, ~20 min for the same counts). Lesson: read my own notebook and STATE T10 line first; the agenda can lag the ledger.
- Traps: scripts that exec the 087 prelude lose `A_` and `P`; guard `standardise` on non-line quivers; A8 build needs a rerun after exit 2.
- Believe: nothing more is learnt by larger F/R radius for this pair; no cross-component control exists.

## Next
- Extend `lnaMoves.discoverMoves` to take a seed row (2400230) and report why no verified move exists (small, with test).
- 3304000 depth-7 J != 0 scan (about 4 x depth-6 cost) only widens a MISS; do it only if asked.
- Find any non-move-set certificate for a non-LNA pair; that is the thing T10 waits for.
