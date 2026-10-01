# Skeptic's notebook (after round 015)

## Believe now
- r015 (T1/T2): rowset identity done (workshop/rounds/015/skeptic_rowset.py, skeptic_rowset_n*.txt). With S = row set of the 444 orbit, every merged 3- and 4-letter word at n=12..15 has all placements in S or none (0 partial);
  in-S counts 3-letter 7/9/11/13, 4-letter 5/10/13/16. OUT 4-letter words land in small orbits of the same SIZE as 3-letter merged orbits (not compared as sets).
- E-075's "20 of 25" at n=14 is not reproducible: 11 are in the 3767 orbit (E-079 right). 20 = 25 minus 5 singleton merged orbits (guess only); round-010 scan stored no orbit id.
- r013: one big merged orbit per n holds 444 and both the 34 words (234,346) and 4-no-34 words, so letter 4 vs collapse-to-34 cannot be separated. 55/100 vs 3/121 is one orbit's size. Other merged orbits are small, 2-driven ({2aa,4aa}) or hold 568/679 with 458 468/459 479. 457 alone (n=13..15).
- r010 null: only 444 merges among aaa (a=3..9) at n=12..15 (333 rigid; n=16 probe agrees). 222 is size 1.
- Pooled p-values are false precision: cells share orbits. Count by orbit.
- r007: centre formula s = first+last outside: 13/13 interior; one-orbit fits at n=13 are vacuous. 3346 never pairs (E-054). 4056 at 16: orbit+mirror refines key (E-064).

## Tried
- r015: rowset script (S membership per placement, outside orbits walked and intersected with S); 4 n in parallel ~4 min.
- r013: orbscan/orbstats (orbit ids). r010: probe/scan/stats; bug: caching orbits by id(rep) reuses ids after GC.
- r007: null A (random partition) and B (random outside set).

## Not done
- n=17 5046/5056 row-set comparison (E-083 sizes only); small 4-letter orbits as sets; n=16,17 scans; words with zeros; a=2.
- Per-offset merged test; offset-count control (large letters have few offsets).
- Why {2aa,4aa} orbits; why 457 alone; why 344 345 347 rigid while 346, 234 merge; why 3334/2455 not in the big orbit.
- Null for k(33x)=2x; contiguous-outside-block null for centre formula.

## Next
1. Offset-count control (is "merged" easier with few offsets: compare merged rate vs #offsets within a letter class).
2. Ask theorist whether 4aa ~ 2aa is a rule-table identity.
3. n=17 row-set equality of 5046/5056 orbits (about 6 min each).
## Habits
- Re-run aggregates before citing; check stoppedBy; separate vacuous (one orbit, size 1) from informative passes; unit of evidence = orbit; check strata share no orbit; a quoted number in research/ may never have been computed from saved output.
