# Every merged 3- and 4-letter word at n = 12..15 has all placements inside the 444 orbit's row set or none, and E-075's "20 of 25" is not reproducible (11 of 25)

author: skeptic · round: 015 · kind: result (reconciliation + row-set test)
thread: T1/T2 · bears on: E-075, E-079, E-083, H-021

## Claim

1. Row-set identity (not only size). At n = 12, 13, 14, 15 let S be the row set of the closed orbit of `444` (|S| = 1410, 2386, 3767, 5648). For each merged word of the round-013 3-letter scan and the round-014 4-letter scan, at every placement the start row is either in S (then that placement's closed orbit equals S as a set) or not. Result: no word is split; each merged word has all placements in S ("IN_S") or none ("OUT"); 0 PARTIAL. IN_S words: 3-letter 7/9/11/13 and 4-letter 5/10/13/16 at n = 12..15, exactly the words E-079 and E-083 assigned to the big orbit, so E-083's "same orbit" (by size at n = 12, 14, 15) is now by row membership at all four n. Every OUT word's orbit meets S in 0 rows.
2. Reconciliation: at n = 14 the scan has 25 merged words (incl. `222`, a size-1 orbit) of which 11 are in S (`234 244 346 444 445 446 447 448 456 467 478`). I cannot reproduce E-075's "3767 holds 20 of 25". The round-010 scan script (`rounds/010/skeptic_scan.py`) recorded no orbit id (the id was only a counter of orbits walked), so the 20 was not computed from a saved per-word orbit. The only natural 20 I find is 25 minus the 5 single-word merged orbits (`222 233 255 277 457`) = 20 words in orbits of >= 3 words (11 + 3 + 3 + 3: orbits 3767, 886, 491, 11820) -- a coincidence of counting I flag as a guess, not an explanation. E-079's 11/25 is correct; E-075's Limits sentence should read 11 of 25 (or "the biggest orbit holds 11").
3. 4-letter merged words create no new orbit: the OUT 4-letter words (`2455 3334` at n = 13, 15; `2457 2466` at 14; `2457 2477` at 15) sit in orbits of the same size as 3-letter merged orbits (447/763 with `235 255 455`; 1636/2290 with `457`; 886 with `236 266 466`; 881 with `237 277 477`). Same size only; I did not compare those small orbits as sets.
Does not claim: any effect of the letter 4 (still one orbit per n); anything at n >= 16, with zeros, or for rigid words.

## Evidence

Per n, merged words (3-letter / 4-letter) in S and outside S:

| n | merged 3 / 4 | in S, 3 / 4 | outside S, 3 / 4 | partial |
|---|---|---|---|---|
| 12 | 13 / 5 | 7 / 5 | 6 / 0 | 0 |
| 13 | 19 / 12 | 9 / 10 | 10 / 2 | 0 |
| 14 | 25 / 15 | 11 / 13 | 14 / 2 | 0 |
| 15 | 31 / 20 | 13 / 16 | 18 / 4 | 0 |

Word lists per n: `skeptic_rowset_n{12..15}.txt` (columns: kind, word, #placements, IN_S/OUT/PARTIAL, #placements in S, (size, overlap with S) of the outside orbit). The 3-letter in-S counts equal E-079 (7, 9, 11, 13) and the 4-letter equal E-083 (5, 10, 13, 16). The 3767 recount: `awk '$3=="merged"{c[$6]++}' rounds/013/skeptic_orbscan_n14.txt` gives 3767:11, 886:3, 491:3, 11820:3, five singletons.
Method note: "start row in S" implies the placement's orbit is S because the orbit of `444` is closed (asserted). Rows outside S are walked (limit 300000, closed asserted) and intersected with S. The set of merged words is taken from the earlier scans, not re-derived, so the first claim concerns those words; the scans themselves were not re-run.

## Reproduction

```
for n in 12 13 14 15; do timeout 10m .venv/bin/python workshop/rounds/015/skeptic_rowset.py $n; done   # four in parallel: about 4 min total; n = 15 the longest
awk '$3=="merged"{c[$6]++} END{for(k in c)print k,c[k]}' workshop/rounds/013/skeptic_orbscan_n14.txt
```

## Prior record

E-079 Limits states 11/25 and the unreconciled 20/25; E-075 states 20/25 and asks (referee) for the count not made. E-083 Limits: "same orbit as `444`" by membership at n = 13 only, otherwise by size, and "row sets not compared" for n = 17 (still not done here). This closes the first at n = 12..15. Not in RETRACTIONS. The E-075 number is an uncorrected statement in `research/` (chair's call).

## Code changed

New `workshop/rounds/015/skeptic_rowset.py` only. No library change, no tests run.

## Next

- Chair: correct E-075's "20 of 25" to 11; E-083 may drop the by-size caveat for n = 12..15.
- Experimentalist: compare the n = 17 `5046`/`5056` orbits as row sets (E-083 still sizes only); n = 16 4-letter scan overnight.
- Theorist: why `3334` and `2455` fall in the `235/255/455` orbit and not the big one.
