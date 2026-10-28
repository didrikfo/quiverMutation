# Toolsmith notebook

## Rounds 003-090 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s, `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); meet up to sign-aware iso. Hom+Cartan has no NO power, the iso meet is the check.
- T10 n = 9: key-kept perI J != 0 steps exist in both F-010 balls and all rejoin their own class by J = 0 moves (E-234, E-239, E-244). The join97 README misreads its log: 26 353 / 15 401 are expansions, None-key drops 9 / 11 (round 090). Power of the 7 + 7 NO is nil beyond 3304000's move component; the NO says nothing on the J = 0 premise.
- Tools: `087/toolsmith_d6scan.py` (build / parents / scan, `--budget-hours 0.14`, env DEPTH), `090/toolsmith_nonekey.py`.

## Round 093
- Ran 2400230(8) vs 3304000(6) on F/R + lnaMoves jump (`093/toolsmith_unionjoin.py lna|extra|meet`): MISS, no WL collision at all, no common lnaMoves component. Whole job about 20 min (builds 8 + 3 min in parallel).
- Believe: 2400230 is an isolated vertex of the lnaMoves graph (component size 1); 3304000's component has 7 rows; A(8) touches 3 components, B(6) 3, none common. So the E-226 gap sits at 2400230: it has no verified move at all.
- Believe: E-244 index 8427's child is isomorphic to index 8063's child (joined at total 6); "8427 first" is vacuous (iso, not labelled-equal).
- Referee 093: the union MISS is implied by 0 common components + E-238 (not independent); F/R meet at 8+6 is the only new search. A8 build needs the rerun after exit 2 (depth-7 pickle gives 59 LNA, not 75).
- Honest limit: a MISS only for one jump layer, depth 8 / 6, no cross-component control exists (no known pair crosses).
- Traps: scripts that exec the 087 prelude lose `A_` (use another name) and `P` (prime for `rank`); a silent exit means exactly that. Nodes of a ball may not be lines: guard `standardise` (quiver with out/in-degree <= 1, n-1 edges) and try/except.
- Next: `discoverMoves` seeded with 2400230 to find why it has no move; 3304000 depth 7 J != 0 scan; stop repeating radius increases.
