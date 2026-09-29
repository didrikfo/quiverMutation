# Skeptic's notebook (after round 001)

## Believe now
- F-053's exceptions 3346 and 4056 are real features, not artefacts. Orbits closed far below
  cap; Coxeter key (`coxeterTables.lnaCoxeterKey`, gauge-free, walk-free) separates 3346's
  offsets at every n to 40, so it never pairs.
- 4056 is a partial reflection with an onset: pairs sum n-13, last three offsets alone.
- The H-020 failures I checked (350066, 6600066) are the same onset effect at n = 13.
- Coxeter-key partition of offsets == reduced-walk orbit partition in all 13 cases compared
  (45, 3346, 4056, 350066, 6600066 at 13-15). Key gives a proof of separation, only a
  necessary condition for joining.

## Tried
- orbitReport at 13/14/15 for the above; key scans to n = 40; 585-word scan at n = 24.
- Scripts kept in workshop/rounds/001/skeptic_{cmp,cox,scan}.py.

## Not done
- The other four H-020 failures (ledger not in repo; need names from a census run).
- Any mismatch between key class and orbit across the catalogue (the case that would matter).
- The 7 cores with triple key classes (7778, 7789, 7899, 7909, 8078, 8889, 9089): orbits?

## Next
1. Whole catalogue at n = 13, 14: orbit partition vs key partition; report any mismatch.
2. Same for cores where key class has 3 members.
3. H-015: the Coxeter key as the second invariant for the guard (it is already used in verifyMove).
## Habits
- Check `stoppedBy`, not just sizes. Grep research for lnaCoxeterKey before claiming novelty.
- Attribute limits: key = proof of difference only.
