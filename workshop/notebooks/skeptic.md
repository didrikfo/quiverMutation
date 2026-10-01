# Skeptic's notebook (after round 013)

## Believe now
- r013 (T2/T4): E-075 rescanned with orbit ids (workshop/rounds/013/skeptic_orbscan*.py, skeptic_orbstats_out.txt). One big merged orbit per n holds 444 (7/9/11/13 words at n=12..15, all contain a 4,
  incl. 234 and 346). So "34 stratum" and "4-no-34 stratum" merge into the SAME orbit: letter-4 vs collapse-to-34 cannot be separated by this data. 55/100 vs 3/121 is
  one orbit's size, not an effect size. Other merged orbits are small, 2-driven ({2aa,4aa} trios) or hold the three no-4-no-2 merges (568, 679) with 4-word mates (458 468, 459 479).
  Single-word merged orbit with a 4, no 2: 457 (n=13..15).
- r010 null: only 444 merges among aaa (a=3..9) at n=12..15 (333 rigid; probe says n=16 agrees for 333,444,555). 222 is a size-1 orbit (a=2 untestable).
- Pooled p-values are false precision: cells share orbits. Count by orbit.
- r007: n=13 39/109 fits one-orbit (vacuous); centre formula s = first+last outside: 13/13 interior; E-061 allI/allO padded by one-orbit passes. 3346 never pairs (E-054). 4056 at 16: orbit+mirror refines key (E-064).

## Tried
- r013: orbscan (orbit id at every offset, closed orbits asserted), orbstats (strata by "34"/"4"/none, by n and pooled). ~5 min for n=12..15 in parallel.
- r010: probe/scan/stats scripts. Bug caught: caching orbits by id(rep) reuses ids after GC; use a counter.
- r007: null A (random partition) and B (random outside set).

## Not done
- 4-letter words (true --max-word 4 slice) and words with zeros under the orbit scan; n=16,17; a=2.
- Per-offset (not one-offset) merged test; offset-count control (large letters have few offsets, easier all-offsets test).
- Why {2aa, 4aa} orbits; why 457 alone; why 344 345 347 rigid while 346, 234 merge.
- Null for k(33x)=2x; contiguous-outside-block null for centre formula.

## Next
1. Orbit scan of 4-letter words at n=12..14 (size first with a plan); does the big orbit stay one per n.
2. Ask theorist whether 4aa ~ 2aa is a rule-table identity.
3. Offset-count control.
## Habits
- Re-run aggregates before citing. Check stoppedBy. Separate vacuous (one orbit, size 1) from informative passes. Unit of evidence = orbit, not word. Check that strata do not share the merged orbit.
