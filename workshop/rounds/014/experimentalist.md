# Saved n = 17 run: `5046` and `5056` each have two closed orbits (122673 at even offsets, 54266 at odd); at n = 12..15 the 4-letter words with a 4 merge into one big orbit per n, the same one as `444`

author: experimentalist · round: 014 · kind: result (saved reruns; the 4-letter scan is new)
thread: T1/T3, T2/T4 · bears on: E-082, E-079, E-076, E-081, E-077, H-021

## Claim

1. `batch.py orbits 17 --cores 5046,5056` (plan: 2 units, 0 done) closes every orbit under the default limit 1500000. For both words: orbit {0,2,4,6} has 122673 rows, orbit {1,3,5,7} has 54266, no pairs, each holds its own mirror. This reproduces E-082's n = 17 sizes for `5046` and now for `5056` too (E-082's referee had not re-run `5056`), with saved output.
2. 4-letter scan (letters 1..9, nondecreasing, contains a 4, >= 4 placements, orbit of the middle offset, reduced walk; 20/35/56/84 words at n = 12..15; no orbit capped at limit 300000): the merged words (all placements in one orbit) are 5, 12, 15, 20 words. At each n one orbit holds most of them: 5/5, 10/12, 13/15, 16/20 words; the others are small orbits (2455 and 3334 at n = 13, 15; singletons `2457`, `2466`, `2477`).
3. The big orbit is the `444` orbit of E-081: its sizes 1410, 2386, 3767, 5648 at n = 12..15 equal that of `234`/`444` in round 013 (`skeptic_orbscan_n*.txt`) (and 4-letter words 2234 2244 2444 2445 2446 4445 ... are in it). Its 4-letter words always carry a 2 or a 4-run (`2xxx`, `4445`, `4456`, `4467`, `4478`, `4566`), so the 4-letter slice adds members to the same one orbit and no new big orbit.
Does not claim: that this is a rate or an effect of the letter 4 (counted by orbit it is one event per n, as in E-081); anything for words without a 4, with zeros, n >= 16, or `--max-word 5`.

## Evidence

n = 17 (`experimentalist_n17_5046.txt`, `_5056.txt`; the two ran in parallel and share one ledger, so the first summary lists both words):

| word | orbit offsets | size | closed |
|---|---|---|---|
| 5046 | {0,2,4,6} | 122673 | yes |
| 5046 | {1,3,5,7} | 54266 | yes |
| 5056 | {0,2,4,6} | 122673 | yes |
| 5056 | {1,3,5,7} | 54266 | yes |

Wall time 5 m 36 s and 5 m 32 s with two running at once (user 3 m 13 s, 3 m 09 s). Identical row-set sizes for both words match E-082 ("same row sets for both words"); I compared sizes only, not the row sets themselves.

4-letter scan, per n (words / merged words / orbits holding a merged word / size of the biggest merged orbit and its words):

| n | words | merged | merged orbits | big orbit (size, words) | `34`-words merged |
|---|---|---|---|---|---|
| 12 | 20 | 5 | 1 | 1410, 5 | `2234` (1 of 10) |
| 13 | 35 | 12 | 2 | 2386, 10 | `2234 2346 3334` (3 of 15; `3334` is in the 447 orbit with `2455`) |
| 14 | 56 | 15 | 3 | 3767, 13 | `2234 2346` (2 of 21) |
| 15 | 84 | 20 | 4 | 5648, 16 | `2234 2346 3334` (3 of 28) |

Orbit lists with words: `experimentalist_orbstats4_out.txt`; raw rows (word, #placements, merged/rigid, #held, orbit id, size, closed): `experimentalist_orbscan4_n{12..15}.txt`. Merged `34` words: no `34` word with 3 letters of `3`/`4` other than `3334` merges; `2234`, `2346` are in the big orbit (they contain `234`/`346`, the E-081 members).
Observation, not tested: `3334` (merged at 13, 15) is in a small orbit (with `2455`) and not in the big one, so the `34` route and the big orbit are again not the same event.

## Reproduction

```
timeout 10m .venv/bin/python batch.py orbits 17 --cores 5046,5056 --plan      # seconds
timeout 10m .venv/bin/python batch.py orbits 17 --cores 5046                    # 5.5 min (alone, about 3.3 min cpu)
timeout 10m .venv/bin/python batch.py orbits 17 --cores 5056                    # 5.5 min, resumes ledger logs/orbits-n17-w4a6-o1500000.jsonl
for n in 12 13 14 15; do timeout 10m .venv/bin/python workshop/rounds/014/experimentalist_orbscan4.py $n; done   # 9 s, 16 s, 58 s, 3 m 10 s
.venv/bin/python workshop/rounds/014/experimentalist_orbstats4.py > workshop/rounds/014/experimentalist_orbstats4_out.txt
```
(Run the two n = 17 commands in parallel as I did, or one `--cores 5046,5056`; `logs/` is not versioned, the saved outputs are in this folder.)

## Prior record

E-082 states the n = 17 sizes (122 673 / 54 266, one run each, `5056` not re-run, no saved file): confirmed and saved here, so the E-082 limit is lifted for `5056`. E-081/E-077 are the 3-letter scan; the 4-letter slice is not in `research/` (grep `max-word 4` and orbit scan: E-076 and E-066 count key-coarser cores, not merged words). Not in RETRACTIONS.

## Code changed

New scripts `experimentalist_orbscan4.py`, `experimentalist_orbstats4.py` only; no library change, no tests run.

## Next

- Theorist: 2234, 2244, 2444, 2445 are in the big orbit at every n: does the rule table reduce each to `244`/`234`? And why is `3334` in a different small orbit with `2455`?
- Skeptic: confirm that the n = 15 big orbit of 5648 is the same row set as the 3-letter one (I compared sizes only).
- Experimentalist later: n = 16 4-letter scan (about 10 min or more, probably overnight with `--jobs`), 4-letter words without a 4.
