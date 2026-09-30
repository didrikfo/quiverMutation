# Skeptic's notebook (after round 007)

## Believe now
- 3346 never pairs (key, E-054). 4056 at 16: orbit-plus-mirror refines key (E-064), "key = orbit" false.
- 6600066 key sums are n-13. My round-001 aggregate was wrong; rerun aggregates before citing.
- NEW r007 null test (workshop/rounds/007/skeptic_null*.py): at n=13, 39 of 109 fits are one-orbit (vacuous);
  70 informative fits: shuffled-partition null fits with P=0.73; only 7/70 < .05, 32/70 < .2. Weak.
- At n=15/16 (12 chosen cores) the null has power: 10/12 and 9/9 < .05. Cores were chosen after fitting at 13.
- Centre formula s = first+last outside: 13/13 interior survives (chance ~4 of 13, joint 1e-7); end-touching 17/21
  vs 2.6 chance; allI/allO counts in E-061 (45/62, 10/13) are padded by 39 one-orbit free passes (informative 9/26, 7/10).
- Null is crude (uniform labels, ignores orbit structure); say so.

## Tried
- r007: null A (random partition, same block sizes, fit()) and null B (random outside set); 300-1000 trials.
- r002: walk 4056 at 16; mirror test; scan at n=24.

## Not done
- Null for k(33x)=2x and for n-independence of k at random cores.
- Fresh cores at 14-15 not pre-selected to test the interior formula (request to experimentalist).
- Contiguous-block null done only by hand estimate; script it.
- The other four H-020 failures (names not in repo); 4056 at 17.

## Next
1. Script a contiguous-outside-block null; add n=14 data when a sample exists.
2. Null for k(33x) across n=13..17.
3. Re-check E-060/E-061 wording once theorist drops the one-orbit cores.
## Habits
- Re-run my aggregates before citing. Check stoppedBy. Key proves difference only. Always separate vacuous passes (one orbit) from informative ones.
