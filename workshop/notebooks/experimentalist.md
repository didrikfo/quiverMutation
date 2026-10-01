# Experimentalist notebook (rewritten each round)

## What I now believe (after round 011)
- H-017 at n = 9: all 16 K = 4 candidates reach nothing at depth 6 (13 new shards in round 011: 308-567 s wall with 4 in parallel on 4 cores; none timed out). Outputs `workshop/rounds/011/experimentalist_cand*.txt`. Bounded negative (E-069: found only if a member within depth).
- Depth 7 would be about 5.5x depth 6, so 25-50 min per candidate: not a 10-minute shard; overnight (Menu 4).
- From 010: key-coarser cores of the `--max-word 4` catalogue: list A (9 words `35 455 3334 3336 5003 5055 5504 5505 5506`) at n = 12, 14, 16; list B (10 words `36 405 466 3335 5004 5006 5046 5056 5066 5605`) at n = 13, 15; orbit+mirror always finer than key, never incomparable. Output `workshop/rounds/010/experimentalist_keycoarser_out.txt`.
- n = 16: the 20300 of 4056 {1,2}, 348 {2,3}, 349 {1,3} is one orbit plus mirror.
- Earlier: k(34x) = x + 3, d = 0; 346 one orbit; 45x no reflection; 4046 reflection k = 11; 5046/5056 translation at odd n; mirror join 344/348/349 at 15..17.

## What I tried
- Round 011: `workshop/rounds/011/experimentalist_shards.sh I` for 13 indices via `xargs -P 4` (about 40 min wall). The launching shell returns at once; wait with a `for ... sleep 20` loop inside one call (timeout up to 590000).
- Round 010: `toolsmith_orbitclass.py N` at n = 14, 15, 16 (12, 22, 60 min with --jobs 4; ledgers in `logs/`, resume windows `timeout 9m`).

## What I would do next
1. Depth 7 at n = 9 for the 16 candidates as overnight (needs a checkpoint or in-candidate budget; ask toolsmith). Time one candidate alone first if a slot allows.
2. K = 100 (160 candidates) at depth 5 is about 2 h: size with --list and one shard.
3. Do the 7 of E-059 (344 366 4044 4403 4404 4405 4605) share orbits across cores at n = 16? Generalise `same20300`.
4. n = 17 key-coarser lists (overnight), `--max-word 5` at n = 14, 34x at n = 18, 44x (E-068 open).
- Watch: a timeout is not a verdict; parallel timings are inflated; 16 candidates are a K cap, not a proof of coverage.
