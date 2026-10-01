# Counted by orbit, E-075's "letter 4" effect is one orbit per n (the `444` orbit), not 55 independent events; words with `34` are in it too

author: skeptic · round: 013 · kind: negative (sharpens E-075; refutes its 55/100-vs-3/121 as an effect size)
thread: T2/T4 · bears on: E-075, E-071, E-065, E-068, H-021

## Claim

Rerun of E-075's scan (all nondecreasing 3-letter words, letters 1..9, >= 4 interior offsets, n = 12..15, size-1 orbits dropped; 291 word-cells) with the
orbit of every word recorded. "Merged" = all offsets of the word lie in the orbit of its middle offset (E-075's test). Findings:
1. At each n there is exactly ONE big merged orbit, the one holding `444`. It holds 7, 9, 11, 13 merged words at n = 12, 13, 14, 15 (all contain a 4).
   Of the merged 4-words, 7/8, 9/12, 11/16, 13/19 sit in it (`34`-words: 2/2 each n, namely `234`, `346`; the rest are `4`-no-`34` words `244 444 445 446 447 448 449 456 467 478 489`).
2. The `34` stratum is therefore NOT separate: `234` and `346` are in the same orbit as `444`. Counted by orbit, "contains `34`" and "has a 4 but no `34`" are the same
   merged event at the big orbit, so this data cannot separate "letter-4 effect" from "collapse to `34`": the two strata share their merged orbit.
3. The other merged orbits are small (1 to 3 words) and are not 4-effects: `{236,266,466}`, `{235,255,455}`, `{237,277,477}`, `{238,288,488}`, `{239,299,499}` hold a 2 (2 is near-free;
   the `4aa` word rides along with `2aa`); `{458,468,568}`, `{459,479,679}` carry the three "no 4, no 2" merges of E-075 (`568` at n = 14, 15; `679` at 15) with 4-words as orbit-mates;
   `233 255 266 277 288 457` are single-word merged orbits (only `457` has a 4, no 2: n = 13, 14, 15; `233 255..` have a 2).
4. What survives of E-075: the observation (`444` the only merged `aaa`, a = 3..9) and "merged words are almost all in 4-or-2 orbits". What does not: the 55/100 vs 3/121 contrast as a
   rate of a letter property. Orbit-collapsed, the same contrast is, over the 4 n pooled (an orbit counted once per n):

| stratum (nondegenerate) | words merged/total | orbits with a merged word / orbits hosting a word |
|---|---|---|
| contains `34` | 8/26 | 4/20 |
| 4, no `34` | 47/74 | 17/33 |
| no 4, has 2 | 26/70 | 18/50 |
| no 4, no 2 | 3/121 | 3/55 (all three in orbits with 4-words) |

   (Strata overlap in orbits, so orbit counts across rows are not additive; and the orbit unit is still weak: most of the 17 "4, no 34" orbits with a merged word are the same
   `{4aa, 2aa}` and 4-in-big-orbit types, 5 to 6 orbits per n.) The orbit-level rate for "4, no 34" (17/33) is not independent evidence either: 4 of the 17 are the big orbit, counted once per n.
5. Not separated, and not claimable: why `344 345 347 348 349` (n = 14, 15) are rigid while `346`, `234` merge; `34` is not sufficient (as E-075 and the round-010 review say), and it is not
   necessary either (`444` orbit holds 4-no-`34` words). The honest statement is: the merged set at each n is {the `444` orbit} plus small 2-driven orbits; the `34` route to it is neither
   confirmed nor refuted by word rates.

Counts per n (words/merged, orbits hosting/with merged): `34`: 5/2 (4/1), 6/2 (4/1), 7/2 (6/1), 8/2 (6/1) at n = 12..15. 4-no-`34`: 10/6 (5/2), 15/10 (7/4), 21/14 (9/5), 28/17 (12/6). No 4: 19/4 (11/3), 34/6 (17/4), 55/8 (22/6), 83/11 (32/8).

## Evidence

Full table, merged orbits with word lists per n: `workshop/rounds/013/skeptic_orbstats_out.txt`; raw scans `skeptic_orbscan_n12..15.txt` (word, offsets, merged/rigid, held, orbit id (per n), orbit size).
The scan has 295 cells (35/56/84/120 per n, equal to the round-010 files without their "done" lines); dropping the 4 size-1 cells (`222`) leaves 291, and the merged count is 84, as in E-075.
Words with zeros, 4-letter words, n >= 16 are not covered ("`--max-word 4`" cores of length 4 not scanned: this is the 3-letter slice of the catalogue only).

## Reproduction

From the repository root:
```
for n in 12 13 14 15; do timeout 10m .venv/bin/python workshop/rounds/013/skeptic_orbscan.py $n; done   # ~5 min total in parallel, n=15 about 5 min alone
.venv/bin/python workshop/rounds/013/skeptic_orbstats.py > workshop/rounds/013/skeptic_orbstats_out.txt  # seconds
```

## Prior record

E-075 (limits: orbit sharing, count of merged words in the `444` orbit "not made": made here). E-065 (4): `44x` lies in one orbit with `333@0`; E-068: `34x` pairs, `346`. Not in RETRACTIONS. The
orbit-collapsed table and the 2-driven orbit list are new.

## Code changed

New scripts only (`skeptic_orbscan.py`, `skeptic_orbstats.py`); no library change, no tests.

## Next

- Do not write "letter 4 is special" into H-021 as a rate; write "all merged words with no 2 lie in the `444` orbit apart from `457` (single-word orbit) and orbit-mates `459 479 679 458 468 568`".
- Theorist: the 2-driven small orbits (`2aa`, `4aa` ~ `2aa`) are a separate pattern worth a rule (is `4aa` = `2aa` an identity of the rule table?).
- Experimentalist: 4-letter words at n = 12..15 under the same orbit scan (`--max-word 4` proper) to see if the big orbit stays one per n; ask the question per orbit, not word.
