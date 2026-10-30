# T7 / H-010 at pair (8:4)(10:4), k = 3, margin 2: 9 of 9 cells stay at overlap 2, but 0 of 4 informative controls lowers it, so the cells are uninformative about "needs a run of three"

author: experimentalist · round: 099 · kind: negative
thread: T7 · bears on: H-010, F-022, E-230, E-241, E-260
scope: ONE new pair shape (8:4)(10:4) (base overlap 2, shared arrows 10, 11), k = 3, margin 2, 6 two-bystander cells with m = 2 plus 3 with an m = 3 bystander (9 cells); controls single bystander at k = 3 (6 cells) and k = 4 (4 cells), of which 4 are informative; plus 3 start-9 cells at k = 3. Capped sample, no proof.

## Response to referee

1. Counts fixed: 6 + 3 = 9 cells (title, scope, Claim); range is 7 to 20 LNAs; "both" is now "three" m = 3 cells sharing 2 arrows.
2. Title and Claim now say 0 of 4 informative controls ((7,3), (12,3) at k = 3 and k = 4). The 6 non-informative ones: (6,4), (12,4) at k = 3 reach 0 LNAs; (12,2) does not embed; (7,4) base 3 at k = 3; (6,3), (13,3) at k = 4 reach 0 LNAs.
3. Start-9 claim run: `timeout 300 .venv/bin/python workshop/rounds/099/experimentalist_single4.py 3 2 "9,2;9,3;9,4" --jobs 4` (output in `experimentalist_start9_k3.txt`, 30 s): (9,2) and (9,3) do not embed; (9,4) has base overlap 3, reaches 4 LNAs, min 2 (2 lowered, from 3 to 2, never below the pair's own 2). So the sentence stands as: "start 9 gives base 3 or does not embed"; the lowering 3 -> 2 is not a lowering of this pair's overlap.
4. Prior record now cites E-029 and the H-010 isolated-pair restriction (HYPOTHESES.md:1023); the open question is dropped.
5. Stepped chain (starts 8, 10, 12) control: not done this sitting, left for next round (needs a new three-relation setup, not a one-line cell).
Also: the closing "no analogue here" is narrowed to "no analogue found among tested controls; the stepped chain is untested".

## Claim

Hit criterion (stated before running): with pair (8:4)(10:4) in place of (8:3)(9:3), a two-bystander cell (each bystander sharing <= 1 arrow with the pair) is a HIT for H-010 if no reached LNA has maximum overlap below 2, and the shape is only informative if a run-of-three-type control (a bystander sharing >= 2 arrows) lowers it to 0 or 1; a MISS is any reached LNA below 2 without such a bystander.

Result: no MISS. All 6 m = 2 cells and the three cells with an m = 3 bystander sharing 2 arrows reach 7 to 20 LNAs, minimum overlap 2. But the controls fail: none of the 4 informative single-bystander controls, (7,3) and (12,3) at k = 3 and k = 4, lowers the overlap (the other 6 are vacuous, do not embed, or have base 3). So for this pair the E-230/E-241 contrast "bystanders alone stay at 2, a run of three goes to 0" cannot be reproduced; H-010's "run of three" (consecutive starts, as (8:3)(9:3)(10:3)) has no analogue among the tested controls with pair starts 8 and 10 (a start at 9 gives base overlap 3 or does not embed; the stepped chain 8, 10, 12 is untested). The honest reading: the 9 cells are consistent with H-010 but the pair shape is not a test that could have failed at k <= 4, margin 2. Not claimed: that no lowering exists for this pair at larger k.

## Evidence

Two-bystander cells (k = 3, margin 2; base overlap 2 everywhere), `experimentalist_pairs_k3.txt`:

| bystanders | shared arrows with pair | reached LNAs | min overlap | secs |
|---|---|---|---|---|
| (6,2)(7,2) | 0, 1 | 7 | 2 | 47 (jobs 1) |
| (7,2)(13,2) | 1, 1 | 13 | 2 | 47 |
| (14,2)(15,2) | 0, 0 | 10 | 2 | 84 |
| (7,2)(15,2) | 1, 0 | 17 | 2 | 104 |
| (6,2)(14,2) | 0, 0 | 18 | 2 | 120 |
| (6,2)(15,2) | 0, 0 | 20 | 2 | 168 |
| (7,3)(14,2) (shares 2) | 2, 0 | 9 | 2 | 102 |
| (12,3)(5,2) (shares 2) | 2, 0 | 9 | 2 | 131 |
| (12,3)(15,2) (shares 2) | 2, 0 | 10 | 2 | 105 |

Single-bystander controls, pair (8:4)(10:4): k = 3 margin 2: (7,3) reached 2, min 2; (12,3) reached 2, min 2; (7,4) base 3, min 3; (6,4), (12,4) reached 0 LNAs (vacuous); (12,2) does not embed. k = 4 margin 2: (7,3) 6 LNAs min 2 (397 s); (12,3) 6 LNAs min 2 (424 s); (6,3), (13,3) 0 LNAs (491, 516 s). Lowered LNAs: 0 in every row. Compare E-230: pair (8:3)(9:3) with a run-of-three bystander reaches 0 at k = 3.

## Reproduction

```
.venv/bin/python workshop/rounds/099/experimentalist_pairs.py 3 2 "8,4:10,4" "6,2:7,2;7,2:13,2;14,2:15,2;7,2:15,2;6,2:14,2;6,2:15,2;7,3:14,2;12,3:5,2;12,3:15,2" --plan
timeout 10m .venv/bin/python workshop/rounds/099/experimentalist_pairs.py 3 2 "8,4:10,4" "7,2:13,2;14,2:15,2;7,2:15,2;6,2:14,2;6,2:15,2;7,3:14,2;12,3:5,2;12,3:15,2" --jobs 4   # 3 m 56 s wall
timeout 10m .venv/bin/python workshop/rounds/099/experimentalist_pairs.py 3 2 "8,4:10,4" "6,2:7,2" --jobs 1   # 48 s
timeout 10m .venv/bin/python workshop/rounds/099/experimentalist_single.py 3 2 --jobs 4        # about 1 min
timeout 10m .venv/bin/python workshop/rounds/099/experimentalist_single4.py 4 2 "7,3;12,3;6,3;13,3" --jobs 4   # 8.6 min
```

## Prior record

E-152 sampled other pair shapes at k = 3 with one bystander only (`theorist_pairshapes_k3.txt`: (8:3)(9:4), (8:4)(9:4)); E-230, E-241, E-260 are all pair (8:3)(9:3). Two-bystander with a pair whose starts differ by 2 is new (no hit in research/ for 8,4:10,4 or start-gap). H-010 (HYPOTHESES.md:1023) is already restricted to an isolated pair, and 2026-09-17 notes that with bystanders an interior move does lower overlap (E-029, 1416 of 7164 applications at n = 10); the "run of three" wording is E-152's finding, so the >= 2 arrows vs run-of-three question is settled there, not open here. F-022 says a pair lowers only with a run of three; this pair has no run of three among the tested controls, which is the finding.

## Code changed

New scripts only: `experimentalist_pairs.py` (pair as a parameter, wrapper over `lnaMoves.localMutationSequences`), `experimentalist_single.py`, `experimentalist_single4.py`. No library change, no tests run.

## Next

Find a lowering control for a start-gap-2 pair first: three relations with starts 8, 10, 12 (a stepped chain), or a bystander at (9,x) of another length; only then are two-bystander cells for it meaningful. m >= 3 at k = 4 for (8:3)(9:3) still costs more than 10 min per cell; k = 5 margin 2 stays an overnight proposal.
