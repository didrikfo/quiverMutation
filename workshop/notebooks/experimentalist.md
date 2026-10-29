# Experimentalist notebook (rewritten round 097)

## What I now believe
- 097 (S-1 item 3): full n = 15 key scan done (2 674 440 rows, 96 shards of ~2 min, four workers, ~46 min wall). Key class of lone 3 (5,6)/(6,5): 18 328 rows; its free-end images at n = 14: 2 keys for K >= 3, 4, 5; 1 key for K >= 6 (K0 = 5 confirmed at class-key level, split not carried by the lone 3). Class (4,7)/(7,4): 21 045 rows, 12 image keys at K >= 5: a messy class. Key level only; orbit/class level open.
- Two-pass plan beats E-224's: scan for target key first, delete only class members (16 s). No n = 14 table needed.
- Row convention: n-vertex LNA is an n - 2 entry row. 094 (T7): k = 5 margin 1 pair (8:3)(9:3) 5 of 5 m = 2 cells overlap 2; margin 2 at k = 5 does not fit 10 min.
- 091 (T10): 2400230 ball sizes 1, 9, 40, 130, 369, 1001, 2691, 7261, 20210; no bridge to 3304000.
- 089 (S-1): |orbit(B_n)| = n 2^(n-5), n = 6..16. 086: k = 4 margin 2 cells 175-530 s. 083: `derivedOrbits(n)` 82 s n = 13, 319 s n = 14. 078: n = 10 90000000 / 00000030 MISS at 5 + 4.

## What I tried
- 097: `rounds/097/experimentalist_scan15.py` (plan, scan, ends), `_batch.sh` (resumable, 4 workers, 9 min), `_split.py`. Hits archived as tar.gz (140 KB).
- Earlier: `k5`, `k5par`, `jscan`, `orb13`, `k4`, `_size`, `_classlaw`, `_jkey`, `_join97`, `_join10`, `_join9`.
- Gotchas: the machine is shared with other personas (4 jobs got ~0.5 core each, 24 shards of 111k rows did not fit, 96 did); `pkill -f` inside a script kills the harness shell (exit 144); a 10 min timeout kills silently, a missing row is not a verdict; write the result at shard end so restarts are idempotent; `removeVertex` can return None (once in the (4,7) class); no pickles of earlier balls survive.

## What I would do next
1. Move-orbit count of the 18 328 rows (E-146 `ends` union-find, `lnaMoves.rewritesOf`, feasible) and whether key 0 / key 1 ends are one `derivedOrbits(14)` orbit.
2. S-1 item 1 at n = 10 -> 9 (not counted; E-194 did 9 -> 8).
3. Overnight: T7 k = 5 margin 2; 2400230 ball depth 9-10; 8427 join.
