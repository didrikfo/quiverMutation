# The H-017 search passes its L = 5 control (42/42 at n = 6; 8/8 at the first 8 near-trivial LNAs at n = 7), and at n = 9 four candidates reach nothing at depth 5

author: maverick · round: 009 · kind: result
thread: T6 · bears on: H-017, E-071, E-065

## Claim

Speculation level: tested on small cases. (1) The verify-style search (`search.linesReachedFrom`, Coxeter guard on) returns the source LNA from a quipu-with-relations member whose recorded walk path has length L = 5: 42 of 42 at n = 6 (every LNA of `ct.lnaStatus(6)`, one member each, longest relation set), 8 of 8 at n = 7 (first 8 LNAs only). At depth L - 1 = 4 it returns the source never (0 of 50 cases). So the E-071 flip at L holds at L = 5. (2) The recorded paths are shortest: for each of the 50 members the source walked at depth 4 does not reach its certificate (50/50; expected, since `reachedQuipuAlgebras` is an exhaustive DFS without `Visited` pruning and keeps the minimum length, so this checks consistency, not an independent BFS). (3) Sizing for n = 9 below-diagonal candidates (see table): about 5.5 to 5.7 times per depth level; depth 5 for the 4 cells of K = 1 takes 5.3 min; depth 6 about 20 min for 4 (one measured, 280 s); depth 7 about 25 min per candidate. It does not claim: no n = 9 negative was obtained at depth 6 beyond the first candidate; 124 of 132 n = 7 LNAs are untested; the control tests inverse-move handling, not independent discovery (E-071's caveat stands); a class with no hereditary member is still uncontrolled.

## Evidence

Control (script builds members from the walk, rebuilds each as a path algebra, runs the search from it):

| n | walk depth / L | LNAs | members (1 per LNA) | hit at L | hit at L-1 | recorded path shortest |
|---|---|---|---|---|---|---|
| 6 | 5 / 5 | 42 (all) | 42 (rels 2:2, 3:31, 4:9) | 42 | 0 | 42 |
| 7 | 5 / 5 | 8 (first 8 of 132) | 8 (rels 3:1, 4:7) | 8 | 0 | 8 |

n = 6 took 175 s, n = 7 (8 LNAs) 111 s, about 7 to 24 s per member. For comparison E-071 had no L >= 5 case. The path-length census printed by `--plan` shows L = 5 members in every LNA at n = 6 (e.g. `0400`: 5, 9, 20, 21, 36 members at L = 1..5).

Sizing at n = 9, `maverick_verify.py 9 D 1 -1` (4 cells, K = 1, the below-diagonal Smith-profile survivors; wall-clock on a shared machine):

| depth | cand 1 (3,1) | cand 2 (3,2) | cand 3 (3,2) | cand 4 (2,1) | reached |
|---|---|---|---|---|---|
| 4 | 9 s | 14 s | 16 s | 14 s | none |
| 5 | 49 s | 79 s | 98 s | 86 s | none (new: depth 5 negative for these 4) |
| 6 | 280 s | not run (cap) | | | none for cand 1 |

Growth 5.4 to 5.7 per level; extrapolation depth 6 about 7 min for the worst candidate (16 candidates: about 1.5 h, in shards of 2 per 10-minute command), depth 7 about 27 min per candidate (16 candidates about 7 h). The forward walks (E-065) needed depth 6 for `3033030`, `4444400`, with path length up to 6, so a depth-6 search from a candidate is the smallest one with resolution matching the forward walk's; depth 5 alone matches only the other seven LNAs.

Side check: forward walk `maverick_reached.py 9 5 4444400` took 42 s and gives the recorded ten-pair set again (no pair on or below the diagonal).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/009/maverick_control5.py 6 5 5 1      # 42/42, 175 s
timeout 10m .venv/bin/python workshop/rounds/009/maverick_control5.py 7 5 5 1 8    # 8/8, 111 s
timeout 10m .venv/bin/python workshop/rounds/009/maverick_control5.py 6 5 5 2 --plan   # path-length census, 89 s
timeout 10m .venv/bin/python workshop/rounds/004/maverick_verify.py 9 5 1 -1       # 321 s, reached [] x4
timeout 10m .venv/bin/python workshop/rounds/004/maverick_verify.py 9 6 1 -1       # cand 1 done in 280 s, hits cap at cand 2-3
```
Outputs: `workshop/rounds/009/maverick_control5_n6.txt`, `maverick_control5_n7.txt`.

## Prior record

E-071 (L = 4 at n = 7, L = 3 at n = 6; "no L >= 5, L shortest not checked": both now filled for n = 6, partly for n = 7); E-065 (forward walks depth 4..6). The depth-4 candidate negative is E-065's "not evidence"; the depth-5 negative for four cells is new but, with the control now passing at L = 5, it excludes only members within 5 steps.

## Code changed

New `workshop/rounds/009/maverick_control5.py` only (extension of round 007's control: selects path length L, adds `--plan`, a shortest check, a depth L-1 negative). No library change, no tests run.

## Next

Overnight is not forced: depth 6 for all 16 candidates is about 1.5 h, i.e. 10 shards of 10-minute commands at 2 candidates each, which the toolsmith could run in a chair slot with `--budget-hours`; only depth 7 (about 7 h) needs `OVERNIGHT.md`, and E-065 already has the n = 9 depth-7 run approved in Menu 4. Toolsmith: add a candidate-index argument to `maverick_verify.py` so shards can split the 16. Experimentalist: run the n = 7 control for the other 124 LNAs (about 30 min, shardable) if a full L = 5 table is wanted.


## Chair note (round 009, after referee)
Sizing corrected: depth 6 for the 4 measured cells is about 30 min (cand 1 measured, 280 s; cands 2-4 scaled 5.7x from depth 5, about 450-560 s each), 16 candidates about 2 h, and one candidate per 10-minute shard (16 commands), not two. Growth per level is 5.4-6.1 (depth 4 to 5), 5.7 from one depth 5-to-6 pair. Claim (2), "recorded paths are shortest", is a consistency check of the same DFS, not evidence. The n = 7 sample is the first 8 LNAs by sort order (near-trivial). The n = 9 depth 5 and 6 outputs were not saved; the 4 measured cells are the K = 1 below-diagonal survivors of 16 candidates, the other 12 were not sized.
