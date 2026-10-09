# The 10 F-041 n = 8 merges have shortest witness length 6 (3 + 3), not <= 4; all 60 edges pass `tiltingPlus` (n = 8, 10 single two-arrow deletions, guarded search to depth 3 per side); `merges.py --witness` now stores paths

author: toolsmith · round: 047 · kind: tool
thread: T10 (iii) · bears on: F-041, H-015, E-148
scope: n = 8; the 10 pairs left open by `lm.ALL_MOVES` (same list as E-148); `search.meetingPoints(..., 3, alsoDual=True)`, i.e. key-guarded gate walks of depth <= 3 per side; one witness per pair (every pair had exactly one shortest meeting). Not covered: n = 9, 10 merges; merges found by `merges.py` itself (no run of it produced a link in my tests); paths outside the guarded graph.

## Claim

(1) `merges.py --witness` (opt-in; default records and return values unchanged) stores in each checkpoint record `witnesses`: for every LNA reached in a linked orbit, the shortest mutation sequence found (vertex labels) and `start` (0 = the member, 1 = its relation dual).
(2) For the 10 F-041 n = 8 merges the shortest witness found has total length 6, split 3 + 3 (meeting at depth 3 from each side); so "depth <= 4" holds per side (3) but not for the whole A-to-B path (6). Since the search covers every split of a total <= 5 into sides <= 3, no witness of length <= 5 exists inside the key-guarded gate graph with duals. It does not say no tilting path of length <= 5 exists outside that graph.
(3) All 60 edges (6 per pair) are gate-admitted, `tiltingPlus` (J = 0) and key-keeping. Direction caveat: the B-half is replayed forward from B, as the search recorded it; the edge actually traversed in the A-to-B route is its inverse, and J = 0 of the inverse step was not tested (E-147: reverse steps are not forward steps of the same algebra).
Refuted by: a pair whose meeting has total < 6, or a recorded half-path whose replay fails the checks (the script asserts replay ends on the meeting key; 10/10 do).

## Evidence

| pair (A -> B) | side A path | side B path | start | edges gate/J0/key |
|---|---|---|---|---|
| 230300 -> 030300 | 1,1,2 | 1,2,1 | member | 6/6/6 |
| 230302 -> 030302 | 1,1,2 | 1,2,1 | member | 6/6/6 |
| 230302 -> 230300 | 8,8,7 | 8,7,8 | dual | 6/6/6 |
| 230400 -> 030400 | 1,1,2 | 1,2,1 | member | 6/6/6 |
| 240030 -> 040030 | 1,1,2 | 1,2,1 | member | 6/6/6 |
| 250002 -> 050002 | 1,1,2 | 1,2,1 | member | 6/6/6 |
| 250002 -> 250000 | 8,8,7 | 8,7,8 | dual | 6/6/6 |
| 304002 -> 304000 | 8,8,7 | 8,7,8 | dual | 6/6/6 |
| 400302 -> 400300 | 8,8,7 | 8,7,8 | dual | 6/6/6 |
| 030302 -> 030300 | 8,8,7 | 8,7,8 | dual | 6/6/6 |

Full output: `workshop/rounds/047/toolsmith_witness.txt`. All ten share one of two path shapes (the deleted relation is at the left end or the right end; the dual runs mirror it), which agrees with the E-148 reading that these merges are one local move, not ten separate facts. The search stores one path per quiver key, so "shortest" is among the first paths found at that length; a tie between different paths of equal length with different edge J-values was not enumerated.
Unit test (`tests/test_merges_witness.py`, n = 5, depth 3): witness off returns `None` and the same reached rows and walk counts as witness on; every stored path replays (gate checked at each step) to the row it names. 27 passed with `tests/test_merge_decisions.py`. An end-to-end `merges.py 9 --depths 3 --witness` run produced no link (expected at depth 3), so the linked-orbit filtering of `witnesses` in `main` is untested by data; `merges.py 6` has no leftover orbits.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/047/toolsmith_witness.py --all      # about 20 s
timeout 10m .venv/bin/python -m pytest -q tests/test_merges_witness.py tests/test_merge_decisions.py -m "not slow"
timeout 10m .venv/bin/python merges.py 10 --depths 5 --witness   # opt-in; adds 'witnesses' to the JSONL
```

## Prior record

E-148 (same 10 pairs: meeting at total 6, all 2396 edges J = 0) and F-041 (meeting depth 3 + 3). New: the actual paths, the replay with per-edge checks along the witness itself (60 edges, not the BFS ball), the total length 6 stated against the "<= 4" question, and the stored-witness option. The result adds little beyond E-148 mathematically (expected: gate = J = 0 on these edges).

## Code changed

`merges.py`: `--witness` flag; `searchFrom` takes an optional 6th element, returns a 7th (`None` unless asked); the record gains `witnesses` only when asked. `tests/test_merges_witness.py` new. Ran that file and `tests/test_merge_decisions.py`. No default behaviour changed.

## Next

Experimentalist/human: deep merges (distance >= 5) with `merges.py <n> --depths 5 6 --witness`, then replay with `toolsmith_witness.py`'s `replay` for per-edge `tiltingPlus` (overnight candidate, not sized). Skeptic: test J = 0 of the inverse steps on the B halves.
