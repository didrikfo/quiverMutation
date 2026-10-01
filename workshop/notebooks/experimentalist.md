# Experimentalist notebook (rewritten each round)

## What I now believe (after round 014)
- n = 17: `5046`, `5056` each have orbits {0,2,4,6} = 122673 and {1,3,5,7} = 54266 (both closed, own mirrors). Saved in `workshop/rounds/014/experimentalist_n17_50{46,56}.txt`. 5.5 min each with two in parallel; `--plan` says 2 units.
- 4-letter words with a 4, orbit scan n = 12..15 (`experimentalist_orbscan4.py`): merged words 5, 12, 15, 20; the big orbit holds 5/5, 10/12, 13/15, 16/20 and is the `444` orbit (size 5648 at n = 15 equals the 3-letter one; only sizes compared). `3334` merges in a small orbit with `2455` at 13, 15.
- From 011: H-017 at n = 9: 16 K = 4 candidates, nothing at depth 6 (bounded negative, E-069). Depth 7 is overnight.
- From 010: key-coarser lists A (even n, 9 words) and B (odd n, 10 words); orbit+mirror always finer than key. k(34x) = x + 3; 346 one orbit; 4046 reflection k = 11.

## What I tried
- Round 014: `batch.py orbits 17 --cores W` twice in parallel (shared ledger: the first summary lists both words). Scan script cloned from skeptic_orbscan.py with a limit arg and a closed flag.
- Round 011: shards by `xargs -P 4`; wait with a sleep loop inside one call (timeout up to 590000).
- Round 010: `toolsmith_orbitclass.py N` at n = 14..16.

## What I would do next
1. Compare row sets (not sizes) of the n = 15 big orbit across the 3- and 4-letter scans; same for `5046` vs `5056` at n = 17.
2. 4-letter scan at n = 16 (time n = 15 was 3 min, so n = 16 about 10 min: use a shard by word range) and 4-letter words with no 4.
3. Depth 7 at n = 9 for the 16 candidates as overnight (needs checkpoint); K = 100 at depth 5 (about 2 h).
4. n = 17 key-coarser lists (overnight), `--max-word 5` at n = 14, 34x at n = 18, 44x (E-068 open).
- Watch: a timeout is not a verdict; parallel timings are inflated; size equality is not set equality; word-level merged counts are not independent events (count by orbit).
