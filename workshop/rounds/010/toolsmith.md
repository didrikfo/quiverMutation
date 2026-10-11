# `toolsmith_verify.py` shards the n = 9 H-017 candidate search; E-074's depth-5 negative reproduces and one candidate at depth 6 fits one 10-minute shard (434 s)

author: toolsmith · round: 010 · kind: tool
thread: T6 · bears on: H-017, E-074, E-065, E-071

## Claim

`workshop/rounds/010/toolsmith_verify.py` is `workshop/rounds/004/maverick_verify.py` (round 009 had no copy of its own; E-074 cites the 004 file) with three options: `--list` (print the numbered candidates, no search), `--cand I[,J..]` (search only those indices) and `--budget-hours H` (exit 2 before starting the next candidate once H hours are spent, naming the candidates not searched). Without options its output is the old script's, plus a `cand I` prefix. It reproduces E-074: the four K = 1 below-diagonal candidates at n = 9 reach nothing at depth 5. One candidate at depth 6 (the K = 1 candidate 2, the one that took longest at depth 5) reaches nothing in 434 s.

It does not claim: that the other 15 (or 159) candidates reach nothing at depth 6; that a budget interrupts a running candidate (it is checked between candidates only, so a shard must be sized so that one candidate fits); anything about classes with no hereditary member.

## Evidence

Candidate numbering depends on K, because K caps candidates per (polynomial, cords, rels) cell and the survivors are enumerated in the same order. K = 1 gives 4 candidates, K = 4 gives 16 (the "16 candidates" of round 009), K = 100 gives 160 (indices 0..159, all survivors with `maxdiag` -1). Indices are stable for a fixed K and n (the enumeration is deterministic; two `--list` runs agreed). So a shard command must name K as well as the index.

| run | depth | candidate (K = 1 index) | result | time |
|---|---|---|---|---|
| 1 | 5 | 0 (cords 3, rels 1) | reached [] | 42 s |
| 1 | 5 | 1 (3, 2) | reached [] | 66 s |
| 1 | 5 | 2 (3, 2) | reached [] | 79 s |
| 1 | 5 | 3 (2, 1) | reached [] | 68 s |
| 2 | 6 | 2 | reached [] | 434 s |

Run 1 totals 4 m 21 s (E-074: 321 s, 49/79/98/86 s; this machine was faster, same answers). Depth 6 / depth 5 for candidate 2 is 5.5, consistent with E-074's 5.4 to 5.7. Applying 5.5 to the other three depth-5 times gives depth 6 of about 230 s, 360 s, 375 s, so all four fit a shard; with the 12 untimed K = 4 candidates the worst case is unknown, and a candidate near the 10-minute cap would be killed by `timeout` with no output. Budget test: `9 4 1 -1 --budget-hours 0.002` ran candidate 0 (7 s), printed `BUDGET SPENT ... not searched: [1, 2, 3]`, exit 2. `--cand 2` ran only that candidate. Depth 6 for the four K = 1 candidates in one command would be about 20 min, over the cap: use one `--cand` per command.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/010/toolsmith_verify.py 9 5 4 -1 --list            # 7 s, 16 candidates
timeout 10m .venv/bin/python workshop/rounds/010/toolsmith_verify.py 9 5 1 -1 --budget-hours 1  # 261 s, reached [] x4, rc 0
timeout 10m .venv/bin/python workshop/rounds/010/toolsmith_verify.py 9 6 1 -1 --cand 2          # 440 s, reached [] ; rc 0
timeout 10m .venv/bin/python workshop/rounds/010/toolsmith_verify.py 9 4 1 -1 --budget-hours 0.002   # rc 2
```
Outputs: `toolsmith_verify_n9_d5.txt`, `toolsmith_verify_n9_d6_cand2.txt`. Shards for all 16 at depth 6: `... 9 6 4 -1 --cand I` for I = 0..15 (about 16 commands, about 2 h by E-074's extrapolation; not yet run).

## Prior record

E-074 (sizing, depth-5 negative, cap at depth 6), E-071 (flip at path length L), E-065. The depth-6 negative for one candidate is not new in kind (E-074 had candidate 1 at depth 6: 280 s); this is a second one with its output saved. The run is a reproduction, not a discovery.

## Code changed

New `workshop/rounds/010/toolsmith_verify.py` only (a copy of the 004 script plus the options; it also puts the working directory on `sys.path` so `families` imports when run from the repository root as a script). No library change; no tests touched.

## Next

Chair/experimentalist: run the 16 depth-6 shards of K = 4 in chair slots (each `timeout 10m`; candidates that die at the cap should be noted and moved to `OVERNIGHT.md`, with depth 7 which is about 27 min each). A cheap improvement if wanted: have `--budget-hours` also stop inside `fm.verify` (needs a hook in `search.linesReachedFrom`; not done because no experiment needs it yet).
