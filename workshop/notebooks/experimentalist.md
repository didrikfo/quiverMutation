# Experimentalist notebook (rewritten each round)

## What I now believe (after round 004)
- The 12 cores of E-060 keep `k = n - s` and signed `d` at n = 13, 14, 15 (12/12) and at 16 for 9/12. The slide rule `s = first o + last o` at 15: 8/12, failures only the known `4045 4506 4556` and all-inside `334`. No new failure of E-060.
- At 16, `3355 3445 46` have no fit: centre still `n - k`, but the middle pair (offsets summing to s) is two separate singleton orbits of equal size 20300 (also seen for `4056` at 16, EXPERIMENTS 2026-09-30). Same signature as the 7 of E-059 at odd n, so the "defect" is an unmerged middle pair, not literally odd n. Parity is core-dependent: these 3 pair at 13..15, fail at 16; the 7 fail at odd n. Not understood; suspect s or n - hi parity or a threshold in n.
- Round 003 results stand (7 cores: pair at even n, strict mirror at odd n, 12..18).

## What I tried
- Ran 24 jobs at once on 4 cores: n = 15 finished slowly (300-570 s each), n = 16 timed out. Machine has 4 cores; run at most 4 census processes. n = 16 census is 160 s per core alone, 12 cores about 12 min at 4 procs; n = 15 is about 40 s alone.
- Data: `workshop/rounds/004/experimentalist_census_12cores_n15.jsonl`, `_n16.jsonl`, slides n15, `experimentalist_shift.py` table.

## What I would do next
1. Cheap: the 3 failures at n = 17 and 18 (unsized; 17 maybe 10 min each). If the census script had `--budget-hours` and `--plan` these would go in OVERNIGHT.
2. Test the unmerged-middle-pair hypothesis on the 7: is the odd-n failure the same thing (centre odd/even)? Tabulate parity of s = n - k against merged/unmerged middle for all 19 cores at n = 12..16.
3. Full 139-core census at n = 12 or 15 still not run (open from round 003).
4. d(c) vs |head - tail| of H-020 is answered for interior blocks (E-060), still open for all-inside cores.
- Watch: a "no fit" here is an unmerged pair with equal size, not a different centre; do not count it as loss of structure without saying so. Census at 16 is near my 10 min limit per core if run in parallel; use xargs -P4.
