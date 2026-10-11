# The n = 7 BFS closure of the both-die classes does not fit in one command and is not showing signs of closing; no hit in 31 000 expansions

author: toolsmith · round: 035 · kind: negative (sizing) + proposal
thread: T5 · bears on: E-126, E-122, H-015

## Claim

The key-preserving BFS of the two n = 7 LNA key classes that carry the 44 both-die candidates (class idx 1: 12 LNAs, 38 targets; class idx 3: 58 LNAs, 6 targets) cannot be closed in one 10-minute command, and a plan run shows no deceleration of the frontier: the per-level growth ratio is 2.4 to 2.5 on class 1 at levels 6 to 8 (not falling) and 2.1 to 1.9 on class 3 (falling slowly). So no shard was run to a verdict. In 240 s per class (11 938 and 19 477 expansions; 29 318 and 40 043 algebras seen) 0 of the 44 targets were reached from an LNA, which repeats E-126(3) at no greater depth (E-126 got 8 and 9 levels in 270 s; the plan run got 8 and 7). It does not claim the classes are infinite (the orbit is finite for fixed n if the key-class mutation graph is, which I did not check), nor that the targets are unreachable: "not found" only.

## Evidence

Script `toolsmith_closure.py plan CLASS SECONDS` runs the same BFS as E-126's `reach` (expand algebras with acyclic quiver; mutate at every gate-admitted vertex; drop children with an illegal relation; keep children whose Coxeter key equals the class key; dedupe by `canonicalKey`) and prints per-level counts. Targets are rebuilt with E-126's `gen` (imported by exec, so the 44 are the same: 38 + 6 for these classes).

Class idx 1 (240 s, 4 cores busy with two jobs):

| level | expanded | s | exp/s | seen | next frontier | ratio |
|---|---|---|---|---|---|---|
| 1 | 12 | 0 | 89 | 56 | 44 | |
| 2 | 44 | 1 | 72 | 154 | 98 | 2.23 |
| 3 | 98 | 1 | 70 | 371 | 217 | 2.21 |
| 4 | 217 | 3 | 65 | 857 | 486 | 2.24 |
| 5 | 486 | 8 | 63 | 1 988 | 1 131 | 2.33 |
| 6 | 1 131 | 21 | 55 | 4 855 | 2 867 | 2.53 |
| 7 | 2 867 | 58 | 49 | 12 033 | 7 178 | 2.50 |
| 8 | 7 083 of 7 178 | 148 | 48 | 29 318 | 17 380 | 2.42 |

Class idx 3: frontier 194, 502, 1 134, 2 434, 5 176, 10 674, 20 566 (ratios 2.59, 2.26, 2.15, 2.13, 2.06, 1.93); 6 levels done by 100 s, level 7 expanded 9 979 of 10 674 at 82 exp/s; seen 40 043.

Sizing. Time per level is frontier / rate, rate 50 (class 1) to 80 (class 3) per second and slowly falling. Class 1, assuming the ratio stays 2.4: level 9 about 6 min, 10 about 14, 11 about 33, 12 about 80 min, 14 about 8 h, and the cumulative seen count passes 10^6 near level 13. Class 3 at ratio 1.9 falling: level 9 about 5 min, 12 about 25 min. There is no level at which the frontier has turned over, so I cannot give a closure time; the honest bound is "more than 10 h for class 1, unknown". One command (10 min) reaches level 9 on class 1, level 8 on class 3 at best, which E-126 already did. Sharding the BFS by first mutation does not cut the cost (the seen set is shared), so `--plan` of the whole closure is the only meaningful shard.

Artefact check (referee's habit): all counts here are the hit test `targets & seen` on `canonicalKey`; the only filters are the key equality and the acyclic-quiver expansion rule, both inherited from E-126 and not changed. If a target needs a path that passes through an algebra with a cyclic quiver, or through a child with a different key (impossible: key is derived-invariant, but `_coxeterKeyOrNone` can return None), the BFS misses it by construction; I did not test how often None occurs.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/035/toolsmith_closure.py plan 1 240   # 4 min; CAPPED, exit 2
timeout 10m .venv/bin/python workshop/rounds/035/toolsmith_closure.py plan 3 240   # 4 min; CAPPED, exit 2
```
Setup (class tables, 44 targets) takes about 30 s inside each run. `run CLASS --budget-hours H` is the same loop with a budget in hours; exits 0 if closed, 2 if the budget is spent.

## Prior record

E-126 (3) reports the capped BFS, 0 of 44, and proposes the closure; STATE round 033/034 lines ask for `--plan` first. No closure of any n = 7 key class is recorded in `research/` (grep "closed", "closure" with n = 7 finds only the E-113/E-126 caps). This is not in RETRACTIONS. The new content is the growth profile, which says the closure is a multi-hour to unknown job, not an afternoon.

## Code changed

New `workshop/rounds/035/toolsmith_closure.py` (plan/run). No library change, no tests touched.

## Next

- OVERNIGHT proposal (Menu: T5): `.venv/bin/python workshop/rounds/035/toolsmith_closure.py run 1 --budget-hours 8` in the background (class 1, 38 targets), second core for class 3 `--budget-hours 4`. Memory: class 1 holds each seen algebra in the frontier; at 10^6 algebras this is about 1-2 GB, check before launching. Honest expectation: it may not close, and then the result is again "capped, 0 hits" with a deeper cap.
- Cheaper alternative, to be sized next: E-126's reverse direction. Mutate each target at all gate-admitted vertices, keep children whose key is unchanged, and look for a child with a linear-Nakayama quiver and relations; a depth-first with a priority for fewer arrows or smaller sum of dim e_iAe_j may reach an LNA without enumerating the class. Toolsmith, next round.
- Experimentalist: the same plan run for n = 7 class 3's six targets shows how a hit would look if the ratio falls below 1.5; it is the class most likely to close.
