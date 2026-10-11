# An n = 8 control with a non-hereditary start (1 to 4 relations) finds its source LNA in 16 of 16 at its own depth and 0 of 13 one step short, at 7e3 to 4e4 nodes, so it is a usable but modest control for the n = 9 depth-6 negatives

author: maverick · round: 014 · kind: tool
thread: T6 · bears on: H-017, E-071, E-074, E-078, E-083

## Claim

At n = 8 (429 LNAs) I built controls whose starting algebra is a quiver with relations (>= 1), not a path algebra: for an LNA, walk out `reachedQuipuAlgebras` to the depth, keep members at recorded path length exactly L with >= 1 relation, rebuild, and search from the member at depth L. The search returns the source LNA in every case: L = 5, 12 of 12 (LNAs 0-5, two members each, all 1 relation); L = 6, 4 of 4 (LNAs 0, 1, 2 with 1 relation; LNA 0 also with 4 relations). One step short (depth L-1) it returns it in 0 of 13 (12 at L = 5, 1 at L = 6). Sizes (undeduped visits / distinct canonical keys):
depth 5, 4.1e3 to 8.9e3 (distinct 0.6e3 to 1.5e3), depth 6, 6.9e3 to 3.85e4 (distinct 0.65e3 to 3.5e3). The n = 9 depth-6 negatives are 5.0e4 to 6.3e4 (E-083). So the control is within a factor of 1.3 to 9 of the negatives (1.3 for the largest control, 38 516 nodes), and the n = 9 searches are the same kind of walk from a non-hereditary start.

It does not claim: that the n = 9 negatives are sufficient (a control tests that the walk undoes a recorded walk, as E-071 said, not discovery); that the controls are typical: all controls with 1 relation have 7 arrows on 8 vertices (a tree), so zero cords, while the n = 9 candidates have 2-3 cords and 1-2 relations; the 4-relation member is also 7 arrows. No control has cords > 0. Only 5 of 429 LNAs were run at L = 5, 3 at L = 6 (10-minute cap; the walk to depth 6 alone is 90-310 s per LNA at n = 8, depending on load).

Answer to the question: usable as a control for the machinery (non-hereditary starts, mutation in both directions, relation-bearing quivers, depth 6, node count of the right order). Not usable as a control for "does depth 6 reach a class with cords"; that needs a member with cords > 0 and relations, which `reachedQuipuAlgebras` members here did not supply at the cheap end. I did not filter for it; a next step.

## Evidence

Script `workshop/rounds/014/maverick_control8.py` (extends 013 `toolsmith_control.py`: MINREL argument, ordering by fewest relations or `HIGH=1` most, `--plan`, `SHORT=1` for depth L-1). Member counts per LNA (L = 6, relations >= 1): LNA 0: 835 (1 rel 214, 2 rels 383, 3 rels 224, 4 rels 14); LNA 2: 423. At L = 5, LNA 0: 677. So non-hereditary members are the majority (plenty to choose from).

| L | LNA | rels | search depth | found | nodes | distinct | s |
|---|---|---|---|---|---|---|---|
| 5 | 0-5, 2 each | 1 | 5 | 12 of 12 | 4 052 to 8 880 (total 72 282) | 624 to 1 512 | 22-41 |
| 5 | same 12 | 1 | 4 | 0 of 12 | 982 to 1 955 (total 16 476) | 265 to 558 | 3-9 |
| 6 | 0 | 1 | 6 | yes | 19 417 | 1 931 | 46 |
| 6 | 1 | 1 | 6 | yes | 26 237 | 2 402 | 61 |
| 6 | 2 | 1 | 6 | yes | 38 516 | 3 534 | 97 |
| 6 | 0 | 4 | 6 | yes | 6 857 | 648 | 22 |
| 6 | 0 | 1 | 5 | no | 4 573 | 815 | 10 |

Timing: search at depth 6 is under 100 s per member at n = 8; the cost of a control is the walk that finds the member (164 s for LNA 0 alone, 91 s for LNA 2), so one control member per LNA is about 4 minutes, and the L = 6 control for all 429 LNAs is about 30 h in shards (not proposed; the control is a sample, and a sample of 5 of 429 is what I have). Node ratios depth 5 to 6 at LNA 0 (same 1-rel member type): 5 080 to 19 417, about 3.8x; the n = 9 depth-6 candidate sizes (50-63k) sit 1.3x to 3x above the n = 8 depth-6 controls; the n = 7 hereditary controls (3.5k to 17k) are lower still. So size grows with n and with the number of relations in the start; a 4-relation start is smaller (6 857), a 1-relation start is larger, since relations restrict mutation.

Reading for the n = 9 negatives: the walk from a relation-bearing start does invert at 2e4 to 4e4 nodes at n = 8, so an n = 9 walk of 5e4 to 6e4 is not suspiciously short or in a regime the walk was never tested in. What remains untested is a start with cords. Also untested: depth 7 (about 2.6e5 nodes at n = 9, E-083; an n = 8 depth-7 control would be about 1e5, 8-10 min of search plus a walk to depth 7 of about 15 min, over the cap).

## Reproduction

All from the repository root.
```
timeout 10m .venv/bin/python workshop/rounds/014/maverick_control8.py 8 5 5 1 2 0 5            # L=5, found 12 of 12, about 5 min plus walks (about 9 min in parallel with other jobs)
SHORT=1 timeout 10m .venv/bin/python workshop/rounds/014/maverick_control8.py 8 5 5 1 2 0 5    # depth 4, found 0 of 12
timeout 10m .venv/bin/python workshop/rounds/014/maverick_control8.py 8 6 6 1 1 0 0            # 3.5 min, found
timeout 10m .venv/bin/python workshop/rounds/014/maverick_control8.py 8 6 6 1 1 1 2            # about 8 min, found 2 of 2
HIGH=1 timeout 10m .venv/bin/python workshop/rounds/014/maverick_control8.py 8 6 6 1 1 0 0     # 4-relation member, found
SHORT=1 timeout 10m .venv/bin/python workshop/rounds/014/maverick_control8.py 8 6 6 1 1 0 0    # depth 5, not found
timeout 10m .venv/bin/python workshop/rounds/014/maverick_control8.py 8 5 5 1 1 0 2 --plan     # member counts only, 85 s
```
Outputs: `maverick_control8_L5.txt`, `_L5_short.txt`, `_L6.txt`, `_L6_short.txt`, `_L6_b.txt`, `_L6_high.txt` in `workshop/rounds/014/` (`_L6_plan.txt` is a mistaken run with L = 1; ignore it).

## Prior record

E-083 / round-013 toolsmith (n = 7 depth-6 control, hereditary-leaning members, node counts, "open: a control with a non-hereditary source at n = 8"): this fills it, at n = 8. E-071 (found iff within depth, n = 7), E-074 (L = 5), E-078 (the negatives). grep of `research/` for "non-hereditary" control at n = 8 found no entry. Consistent with all of them; nothing contradicted.

## Code changed

New `workshop/rounds/014/maverick_control8.py` only; no library change, no pytest run needed. Not cross-checked against a second implementation: the member is rebuilt from the certificate by `algOf` (as in 013).

## Next

Toolsmith: filter control members for cords > 0 (arrows >= n) and relations >= 1 at n = 8, L = 6, so the control matches the n = 9 candidates (2-3 cords); `reachedQuipuAlgebras` members have arrows in the certificate, so it is one condition in the member filter. Chair: the n = 9 depth-6 negatives can be cited as "run at a node count consistent with relation-bearing controls that were found"; depth 7 overnight stays proposed, with an n = 8 depth-7 control of 1-2 members first (about 25 min each, shard).
