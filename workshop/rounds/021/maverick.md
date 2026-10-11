# At n = 8, an LNA has cord members within 3 mutation steps iff one of its relations has at least 3 arrows (365 of 429, 0 mismatches); a non-MONO L = 5 walk of the rest is about 4.5 CPU-hours

author: maverick · round: 021 · kind: result
thread: T6 · bears on: H-017, E-089, E-091, E-094

## Claim

Cord = quiver with arrows >= n and no parallel arrows, in the search of E-089/E-094 (both directions). Criterion C: the n = 8 LNA (relation lengths d_1..d_6, d_i = arrows of the relation starting at vertex i) has a cord member within path length 3 iff max d_i >= 3, i.e. it is not radical-square-zero (all relations of 2 arrows or none). Tested on all 429 n = 8 LNAs at L = 3: 365 predicted yes, 365 have cords; 64 predicted no, 0 have cords. Of the 64 "no", all 64 also have none at L = 5 . No monomial cord at L = 3 for any of the 429 (MONO = 1). Speculation level: tested on small cases (n = 8, L <= 5 for the negatives; L = 3 for the positives). It does NOT claim that the 64 have no cords at any depth, nor the n = 9 statement, nor anything about monomial cords beyond depth 3 (E-094's n = 8 depth 5 MONO run covers six LNAs only).

Heuristic reason (idea): a relation of >= 3 arrows, a path a -> b -> c -> d = 0, is the case where a mutation can replace a monomial by a sum of two paths (the commutativity square of E-094); with only length-2 relations (rad^2 = 0 on pieces) the mutated algebras stay without a pair of parallel paths. Not proved.

## Evidence

Cord members within L = 3 per LNA (counts, walk about 1.5 s each); minimum path length of a cord member over the 365 positives: 1 for 347, 2 for 17, 3 for 1.

| class | LNAs | with cord member (L=3) | with cord member (L=5) |
|---|---|---|---|
| max d_i >= 3 | 365 | 365 | (6 walked in E-094: 6 of 6) |
| max d_i <= 2 | 64 | 0 | 0 of 64 |
| MONO = 1, all | 429 | 0 | not run |

Remaining question for the 64: radical-square-zero LNAs include the hereditary A_8 (000000); 64 = 2^6 sequences over {0,2}.

Sizing (`--plan` equivalent; the walk is the plan): L = 5 walk per cord-bearing LNA takes 30-60 s serial (E-089 plan: 6-36 s on the first 14; the 64 non-cord LNAs took 915 CPU-s in total, max 72 s, with 4 parallel processes). For the other 365: about 365 x 45 s = 4.5 CPU-hours, about 70 min on 4 cores in four shards, each shard under the 10-minute limit only if cut into about 9 pieces of about 10 LNAs. The MONO filter is a flag on the same walk, so one walk gives both counts. That is a proposal for `OVERNIGHT.md`, not run here. What it would answer is less than before: criterion C already says which 365 carry cords, so the new content of the walk is only MONO at L = 5 (E-094 showed none for 6 LNAs) and n = 8 members with more than 9 arrows.

## Reproduction

```
for i in 0 1 2 3; do timeout 10m .venv/bin/python workshop/rounds/021/maverick_predict.py 8 3 $((i*108)) $((i*108+107)) > workshop/rounds/021/maverick_predict_L3_s$i.txt; done   # about 3 min per shard, shards in parallel
.venv/bin/python workshop/rounds/021/maverick_criteria.py workshop/rounds/021/maverick_predict_L3_s*.txt          # criterion table
IDX=<indices of the 64 zero LNAs> timeout 10m .venv/bin/python workshop/rounds/021/maverick_predict.py 8 5 0 0     # L = 5 for the 64 (four shards, about 4 min)
MONO=1 timeout 10m .venv/bin/python workshop/rounds/021/maverick_predict.py 8 3 LO HI                              # monomial only
```
Output files: `maverick_predict_L3_s*.txt`, `maverick_predict_L5_zero_s*.txt`, `maverick_predict_MONO_L3_s*.txt` (lines: index, relation lengths, number of cord members, minimum path length, seconds).

## Prior record

E-089 (cord members exist only for LNAs "with a 3 early in the sequence", indices 4, 9-13 among the first 14), E-094 (all cord members carry a sum relation; no monomial cord at n = 4, 5, 8), E-091 (MONO at L = 6, six LNAs). Not recorded anywhere (grep for "radical square" and "max.*>= 3" in `research/`): the full-429 classification and its exact iff at L = 3. It sharpens E-089's observation from a guess about position to a statement about relation length: the six cord LNAs in indices 0-13 are exactly those with a 3 or 4, and 0-8 and 14 are the rest. It does not touch H-017 itself: the cords here are sum-relation cords from LNAs, which is the control side, not the E-078 candidates.

## Code changed

New: `workshop/rounds/021/maverick_predict.py`, `maverick_criteria.py`. No library code touched; no tests needed.

## Next

- Chair: `OVERNIGHT.md` proposal: L = 5, non-MONO plus MONO count on the 365, in ten-LNA shards (about 70 min on 4 cores).
- Theorist: prove C. One direction: a relation of >= 3 arrows is mutated at a vertex into a sum relation; other: for rad^2-zero LNAs, tilting mutation preserves "no two parallel paths of length >= 2" (check against step 7 of the mutation).
- Experimentalist: does C persist at n = 9 (1430 LNAs; L = 3 time per LNA unmeasured)? The one LNA with minimum cord depth 3 is worth a look: is L = 3 itself the edge of C?
