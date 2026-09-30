# Skeptic's notebook (after round 002)

## Believe now
- 3346 never pairs: distinct Coxeter keys at every n to 40; orbits agree at 13-15. Proof of non-equivalence.
- 4056: key classes {o, n-13-o}, last three alone. Walk = key at 13-15.
- NEW n=16, 4056: walk gives {0,3},{1},{2},...; key gives {1,2}. Orbits of 1 and 2 are mirror images
  (mirrorRow of one start is in the other orbit), so derived equivalent by F-026. Key = orbit is false;
  key = orbit-up-to-mirror is untested beyond this case.
- My round-001 aggregate was wrong: at n=24, 585 words = 309 one-centre, 155 with key class >= 3, 121 none
  (not 7 / 123). Output saved in rounds/002/skeptic_scan_n24.txt.
- 6600066 key sums are n-13.

## Tried
- Walk 4056 at 16 (252 s, limit 1.5M, all closed); mirror test on offsets 1, 2; scan at n=24 rerun.
- Scripts: workshop/rounds/002/skeptic_{cmp,scan,mirror16}.py.

## Not done
- Orbit-up-to-mirror vs key across the catalogue at 14-16; 4056 at 17 (too long for one call?).
- The other four H-020 failures (names not in repo).
- Scan at n other than 24 (onset behaviour across n for the 155/121 words).

## Next
1. Test whether every key class is one orbit or one mirror pair of orbits (13-16, catalogue).
2. Check the F-053 census claim "self-dual" at 16 for other cores.
## Habits
- Re-run my own aggregates before citing them. Check stoppedBy. Key proves difference only.
