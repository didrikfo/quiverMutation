# Experimentalist notebook (rewritten each round)

## What I now believe (after round 007)
- `k(34x) = x + 3` (the footprint), `d = 0` (s = hi), for x = 4, 5, 7, 8, 9 at n = 14..17, 20/20 cells, all closed; `346` never fits (one orbit, n = 14..17). Not `2x`. Raw singletons in the middle are equal-size mirror pairs or the centre (`d_eff = 0` everywhere).
- `45x`: no reflection; 456, 457 always one orbit; 455, 458, 459 split by offset parity at some n (455 even n, 458 odd n, 459 even n) and merge at others.
- `4046` is a reflection (k = 11, d = 1, n = 12..16, middle pairs unmerged but equal size). The period-2 translation is `5046`/`5056` at odd n (13, 15); one orbit at even n (12, 14, 16). E-060's claim that `4046` also gives {0,2},{1,3} at 13 did not reproduce (I get {0,2},{1},{3}).
- Earlier (004) results stand: 12 cores of E-060 keep k, d; "unmerged middle pair" is mirror pair (E-064).

## What I tried
- `workshop/rounds/007/experimentalist_kd.py`, `_table.py`, `_core4046.py`. Timing: 34x n = 16 about 8 min total (347 200 s, 348 170 s), n = 17 347 alone 370 s, 348/349 about 150 s each; 4 procs in parallel fit in 10 min except 34x n17 (rerun per x).
- Orbit-only walks; mirror equality inferred from equal orbit size, not joined.

## What I would do next
1. Mirror-join check (toolsmith_orbitclass.py) of the equal-size singleton pairs in 344/348/349 at 15..17: turns "equal size" into "same orbit".
2. 34x at n = 18 for x = 7..9 (OVERNIGHT, ~40 min each); 44x not run.
3. Check `4046@13` offset 3 vs `5046@13` {1,3} (equal size 2116): same orbit?
4. Still open from 004: 3 failures at n = 17, 18; 139-core census n = 12 / 14 (overnight).
- Watch: pair sums from classes that also contain 3 members are not fits (the table script prints k for 458/459 wrongly; ignore). Sizes at n = 17 reach 120k; cap is 300k.
