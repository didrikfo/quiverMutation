# `longSquare` now sees parallel-arrow long squares (4 of 16 out-degree 1 rejects at n = 8 class 1 were missed), and the n = 9 class 0 reject walk has checkpoint/resume with `--budget-hours`

author: toolsmith · round: 026 · kind: tool
thread: T5 · bears on: E-105, E-107, E-108

## Claim

(1) `longSquare` (only ever in workshop scripts, rounds 022/023; not in `quivermutation/`) compared predecessor *vertices*
`q[-3]` on the vertex-list `alg.rels`; two paths through a doubled arrow project to the same vertex list, so it returned False.
The fixed test reads the arrow relations (`procedure.relationsFrom`, which prefers `arrowRels`) and asks for distinct second-to-last
*arrows*. On an n = 8 class 1 walk (420 s, 10 106 expansions) all 16 distinct-step out-degree 1 rows with J != 0 now have a long
square; the old test missed 4 of the 16 steps (4 distinct (parent, v), all parallel-arrow). It does NOT claim "reject iff long
square" (E-107: out-degree >= 2 rejects have none, by construction of the test) nor anything about n = 9 beyond the prefix below.
(2) `toolsmith_rejwalk.py` is the round 023/025 walk plus classifier with `--ckpt`, `--budget-hours`, `--max-exp`; resume gives
the same table as an uninterrupted run (n = 7 class 0, 3000 expansions, 3 slices; diff empty). Overnight n = 9 class 0 is feasible
as a multi-night resumable job, not a one-shot (sizing below).

## Evidence

Unit tests, `tests/test_longsquare.py` (3 pass): case 1 builds 3 -> 2 => 5 -> 4 (2 => 5 doubled) with the relation
`a g0 e = a g1 e`; `alg.rels == [[[3,2,5,4],[3,2,5,4]]]` (the E-108 shape), `longSquareOld(alg,5)` False (fails before),
`longSquare(alg,5)` True. Case 2: an ordinary square and its two non-squares (relation stops at v; monomial) agree old/new.
Case 3: parallel arrows before the penultimate arrow (same second-to-last arrow) is correctly not a square.
Old vs new over n = 6, 7 LNAs and duals plus one step out (2 016 and 8 316 (algebra, v) tests): 0 differences, but also 0 positives,
so that check is vacuous; the real agreement check is the walk table, where the `old` column equals `longsq` on all rows except
the 4 parallel steps. Walk tables (`longsq`, `old`):

| walk | J != 0, out 1, new True, old True | new True, old False | J != 0, out >= 2 |
|---|---|---|---|
| n = 8 c1, 420 s, 10 106 exp | 12 | 4 | 38 (class D-part 1/2) |
| n = 7 c0, 3 000 exp | 26 | 0 | 0 |
| n = 9 c0, 100 s, 1 388 exp | 32 | 0 | 0 |

Distinct (parent, v) at n = 8 c1: 54 = 38 D-part1/2 + 8 (dim J 2) + 4 sq + 4 sq-parallel-only. Counts depend on load and cap (E-108).

Resume check: `toolsmith_rejwalk_resume_check_uninterrupted.txt` (one run, `--max-exp 3000`) vs
`..._3slices.txt` (`--max-exp 1000`, `2000`, `3000` with one `--ckpt`; checkpoint 4.8 MB): identical after dropping time and slice count.

Sizing for n = 9 class 0 (`--plan`: 2 starting algebras, class 0 is the smallest). Measured: 14 expansions/s, about 1 KB of checkpoint per
distinct algebra (3 364 algebras: 3.2 MB). BFS levels (new algebras per level) 8, 22, 56, 132, 328, 851, ratio about 2.5 and not
yet flattening, so level 12 would hold about 2e5; the walk may not close. n = 7 class 0 did not close in 10 min either (level 9:
17 099 next at 321 s). A 7 h slice is about 350 k expansions, about 0.5 GB of checkpoint; one save takes seconds.

## Reproduction

```
.venv/bin/python -m pytest -q tests/test_longsquare.py -m "not slow"                       # 2 s
timeout 10m .venv/bin/python workshop/rounds/026/toolsmith_rejwalk.py 8 --class 1 --budget-sec 420 --show 0   # 7 min; rows above
timeout 10m .venv/bin/python workshop/rounds/026/toolsmith_rejwalk.py 9 --plan             # seconds
# resume check (about 3 min):
R="workshop/rounds/026/toolsmith_rejwalk.py 7 --class 0"
.venv/bin/python $R --max-exp 3000                                    # exit 2
for m in 1000 2000 3000; do .venv/bin/python $R --max-exp $m --ckpt /tmp/rw7.pkl; done   # exit 2 each; tables equal
```

Proposed overnight (NOT run; human to approve), solo, not parallel, rerun the identical command after each exit 2:

```
.venv/bin/python workshop/rounds/026/toolsmith_rejwalk.py 9 --class 0 --budget-hours 7 --ckpt /tmp/n9c0.pkl > /tmp/n9c0_slice$(date +%s).txt
```

Exit 0 = closed, 2 = budget spent (checkpoint written; SIGTERM also checkpoints), other = crash (the last periodic save, every 900 s, is kept).

## Prior record

E-108 (STATE T5): "out-degree 1 no-square rejects are parallel-arrow long squares missed by `longSquare`" -- now fixed and counted (4 steps
here). E-107: n = 9 c0 prefix (8 231 algebras, 500 s) shows no D' rejects; coverage-dependent, my 1 388-expansion prefix is shorter and
adds nothing. Prior checkpointed walk: round 019 `toolsmith_walk.py` (scholar_walk, different counts). Nothing in RETRACTIONS touched.

## Code changed

New: `tests/test_longsquare.py`; `workshop/rounds/026/toolsmith_longsquare.py` (`longSquare`, `longSquareOld`),
`toolsmith_rejwalk.py`, two resume-check text files. No file in `quivermutation/` changed (the function was never there), no edits to
older scripts, which keep the old test. Tests run: `tests/test_longsquare.py` only (3 passed).

## Limits and referee notes

- `longSquare` now needs a `PathAlgebra` with `arrowRels` for parallel quivers; for hand-built ones with parallel arrows and no
  `arrowRels`, `relationsFrom` lifts `rels` and raises on ambiguity -- by design, not caught.
- Distinctness of second-to-last arrows is weaker than of vertices, so new True is a superset of old True; on the walks it adds only the 4.
- The walk's `rej` stores one algebra per distinct (parent, v) with J != 0 in memory and in the checkpoint.
- Resume equality is shown for deterministic `--max-exp` slices; time slices differ only in where they stop.

## Next

experimentalist: approve/run the overnight command; classify `sq-parallel-only` with a minimality check (E-105 caveat: is the relation minimal?).
theorist: whether the n = 9 level growth (ratio 2.5) means class 0 never closes, which would make "coverage by algebra count" the only honest unit.
