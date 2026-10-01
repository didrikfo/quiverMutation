# "Only a = 4 merges" survives a neighbour-aware null, but only as "the letter 4 is special"; it does not single out the collapse to 34

author: skeptic · round: 010 · kind: negative (partial: observation supported, mechanism not)
thread: T2 · bears on: E-071, E-065, E-068, H-020, H-021

## Claim

E-071 / round-009 theorist: seed `aaa` collapses to `34` only for a = 4, so `44x` merges into one orbit and `aax` (a != 4) stays rigid.
Null: scan EVERY nondecreasing 3-letter word xyz (letters 1..9, >= 4 interior offsets) at n = 12..15 (291 word-cells), call it *merged* if the
reduced-walk orbit of the word at one interior offset holds it at all its offsets, and ask whether `444` stands out among words.
Result: (1) the observation holds: among a = 3..9, only `444` is merged, at all four n (12..15); `333 555 666 777 888 999` are rigid in every cell
(22 cells). (2) But merging is a property of the letter 4 across the board: 55/100 nondegenerate words containing a 4 are merged, against 3/121 of
words with neither 4 nor 2 (and 26/70 with a 2 and no 4; a leading 2 acts as a near-free letter, `222` is a size-1 orbit, so a = 2 is not testable here).
So "a = 4 special" is expected once letter 4 is a translator, and is not evidence for the particular route `444 -> 34`. NOT claimed: that the 34
route is wrong; only that this data cannot tell it from "any word with a 4 merges easily".

## Evidence

Merged rate among nondegenerate words (orbit size > 1), n = 12..15 pooled (`skeptic_stats_out.txt`):

| subset | merged / words |
|---|---|
| all | 84/291 = 0.29 |
| contains a 4 | 55/100 = 0.55 |
| no 4 | 29/191 = 0.15 |
| no 4, no 2 | 3/121 = 0.02 (the 3 are `568` at n=14, 15 and `679` at 15) |
| no 4, has a 2 | 26/70 = 0.37 |
| letter c (no 4), c = 3,5,6,7,8,9 | 0.16, 0.11, 0.12, 0.09, 0.10, 0.11 |

The letter-4 rate is the largest of the eight letters (0.55; next is 2 at 0.35, others 0.20-0.27 when 4 is allowed). Among `aaa`, 1 of 7 (a = 3..9) merges,
and it is the one with the highest letter rate: the claim is a special case of a pattern, not a lone anomaly. A word-level pooled null for "exactly the
4-seed merges" is not computed as a p-value: the cells are heavily dependent (many words share an orbit: 3767 at n = 14 holds 20 of the 25 merged words),
so a p-value would be false precision.

Why the 34 mechanism is not singled out: words with a 4 and no `34` merge (`445 446 447 448 456 457 467 468 478`), while words with `34` often do NOT
(`344 345 347 348 349` rigid at n = 14, 15; `346` merges; cf. E-068, `34x` pairs). If "34 reaches the slider" were the whole mechanism, `34x` should merge;
it does not. The slider route holds for the seed `444` (theorist's 7-step path) but is not what makes a word merge in general.

`aaa` at n = 12..15 (merged?, orbit size): 444: yes at 1410, 2386, 3767, 5648 (orbit size grows, same word); 333 rigid (98, 501, 491, 877); 555 rigid
(1766, 4217, 2179, 18416); 666 rigid; 777 rigid (n >= 13); 888 rigid (n >= 14); 999 rigid (n = 15, only 4 offsets). Small-offset cells (`888`, `999`) are weak: 4-5 offsets.
Check on the cell counts: `777` at n = 14 held {2} only of 5 offsets, `888` held {1,2} of 4: pairs, as the rigid chain structure predicts.

Weaknesses of my null: letters/offset counts are confounded (words with larger letters have fewer offsets, so fewer ways to look rigid... and merged requires
all offsets, an easier test for short words); orbits capped at 300000 never hit (all closed); one interior offset per word (the same choice as the probe, not as
the theorist's per-offset walk); n = 12..15 only; a = 2 degenerate; words with zeros and 4-letter words not scanned.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/010/skeptic_scan.py 14 4     # ~1 min; n = 12, 13, 15 likewise (writes skeptic_scan_n14.txt)
.venv/bin/python workshop/rounds/010/skeptic_stats.py                      # seconds; output in skeptic_stats_out.txt
timeout 10m .venv/bin/python workshop/rounds/010/skeptic_probe.py 14 44 333 444 555   # orbit and held offsets of single words
```
Scan times: n = 12 4 s, 13 25 s, 14 58 s, 15 135 s.

## Prior record

E-071 and the round-009 theorist submission state the criterion; the review asked for exactly this null. `grep` of `research/` for "letter 4" /
"rate" of merging finds nothing; the near-universality of merging for 4-words is new but unrecorded as a finding. E-068 (`34x` pairs, `346` one orbit) agrees
with my `344 345 347 348` rigid, `346` merged. Not in `RETRACTIONS.md`.

## Code changed

None in the library; new scripts `skeptic_probe.py`, `skeptic_scan.py`, `skeptic_stats.py` (no tests touched).

## Next

- Keep the E-071 wording as "only `444` merges among `aaa`, n = 12..15" and do NOT write the "collapse to 34" criterion into H-021's text (STEERING q2): the null
  does not discriminate it; propose the competing statement "words containing a 4 merge at 0.55, others at 0.02-0.15" for the theorist to explain (why 4?).
- theorist: why are `344 345 347` rigid while `346` and `446` merge; what is special about 4 as a letter in the rule table.
- experimentalist: `aaa` at n = 16, 17 (a = 3..9) and 4-letter seeds `aaaa`; `22x` needs a different test (2 is degenerate: 222 orbit size 1).
