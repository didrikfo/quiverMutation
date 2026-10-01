# All 16 n = 9 H-017 candidates reach no class at depth 6: the 13 untried K = 4 shards each finish under 10 minutes and find nothing

author: experimentalist · round: 011 · kind: negative
thread: T6 · bears on: H-017, E-072, E-073

## Claim

At n = 9, for each of the 13 K = 4 candidates not yet searched (indices 1, 2, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15), the depth-6 mutation search of `toolsmith_verify.py` reaches no quipu-class member (`reached []`). With E-072/E-073 (indices 0, 4, 8 done), all 16 K = 4 candidates are now negative at depth 6. None timed out. This extends the H-017 negative from depth 5 to depth 6 for the 16 candidates.

It does not claim: anything at depth 7 or above; anything for the 160 candidates of K = 100 (only the 16 of K = 4 were searched); that H-017 is true. By E-069 the search finds a class only if a member lies within its depth, so this is a bounded negative, not a verdict. Times are wall seconds with four shards running at once on 4 cores, so they are inflated against a lone run (compare E-073: 434 s alone for K = 1 candidate 2).

## Evidence

| cand | cords | rels | polynomial cell | result | seconds |
|---|---|---|---|---|---|
| 1 | 3 | 1 | P^(1,1,1)_(1,0,1,1) | reached [] | 308 |
| 2 | 3 | 1 | same | reached [] | 359 |
| 3 | 3 | 1 | same | reached [] | 524 |
| 5 | 3 | 2 | same | reached [] | 453 |
| 6 | 3 | 2 | same | reached [] | 488 |
| 7 | 3 | 2 | same | reached [] | 314 |
| 9 | 3 | 2 | P^(1,1,2)_(1,0,0,1) | reached [] | 524 |
| 10 | 3 | 2 | same | reached [] | 518 |
| 11 | 3 | 2 | same | reached [] | 567 |
| 12 | 2 | 1 | P^(1,3)_(1,1,1) | reached [] | 474 |
| 13 | 2 | 1 | same | reached [] | 322 |
| 14 | 2 | 1 | same | reached [] | 365 |
| 15 | 2 | 1 | same | reached [] | 309 |

Seconds are shell wall time of the shard (rc 0 for all 13); the script's own per-candidate time is a few seconds lower (cand 1: 296 s). Slowest: cand 11, 567 s, only 33 s below the cap. Unknown index 12 of E-073 ("untimed"): 474 s (the notes' "12" read as an index). No traceback in any output.

## Reproduction

From the repository root; one shard per index, I in {1,2,3,5,6,7,9,10,11,12,13,14,15}:
```
workshop/rounds/011/experimentalist_shards.sh I     # 300-570 s each; runs toolsmith_verify.py 9 6 4 -1 --cand I under timeout 10m
```
I ran them with `xargs -P 4` (about 40 min total). Outputs `workshop/rounds/011/experimentalist_cand{I}.txt`, times `workshop/rounds/011/experimentalist_shard_times.txt`.

## Prior record

E-072 (depth 5 negative for all 16 candidates; depth 6 for candidate 1 = index 4, 280 s), E-073 (K = 1 candidate 2; K = 4 indices 0, 4, 8). Not previously recorded: depth 6 for the other 13. Nothing in `research/RETRACTIONS.md` touched (not checked beyond the E-numbers above).

## Code changed

New `workshop/rounds/011/experimentalist_shards.sh` (wrapper). No library change, no tests.

## Next

- Chair: record in EXPERIMENTS as an E-entry; T6 can state "16 of 16 K = 4 candidates negative at depth 6".
- Depth 7 is about 5.5 times depth 6: 25 to 50 min per candidate (alone, probably 25 to 40 min), so a `timeout 10m` shard cannot hold it. Overnight proposal (already Menu 4): 16 candidates x about 30 min on 4 cores is about 2 to 3 h.
- Toolsmith: `--budget-hours` does not stop inside a candidate; a depth-7 run needs that hook or a checkpoint.
- Skeptic: the K = 100 set has 160 candidates; K = 4 is a cap, not a sample proven representative.
