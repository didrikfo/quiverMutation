# `compare2` now asserts dim K Q/I' = dim End(T); every dropped relation is rejected (125 of 125, incl. step 12 and class 2), while the 25 true children still pass

author: toolsmith · round: 059 · kind: tool
thread: T10 · bears on: E-174, E-177
scope: n = 7, classes 1 (16 failing steps) and 2 (9 failing steps), label-preserving; wrong algebras = the child's own crels with ONE relation deleted (125 variants); n = 10 and path edges not rerun; equal-dim wrong algebras not re-done (E-177 covers killed-line variants at 7 c1 steps, not step 12 / class 2).

## Claim

In `workshop/rounds/057/toolsmith_endt2.py`, `compare2` now computes dim K Q/I' per vertex pair from `crels` (`quotient_dims`, linear algebra on all paths; `dimcheck`) and returns `NO (dim K Q/I' != dim End(T))` before the Groebner test if it differs from the child's dimension matrix. At all 25 failing n = 7 steps (c1 16, c2 9) the assertion passes for the true child and the verdict is still `iso` (unchanged from E-174). The control E-177 left open is now settled: deleting any one relation of `crels` raises dim K Q/I' by 1 to 4 in 125 of 125 cases (c1 86 incl. all 7 for step 12, the two-parallel-pair step; c2 39), `dimcheck` rejects 125 of 125, while `symcheck2` alone said `iso` 125 of 125. It does NOT claim: that `crels` contains no wrong *extra* relation of the right dimension (that is E-177's killed-line family, done only for c1 single-pair steps), that the complete-relation-set input is verified independently of `child_info` (the dimension is compared with the child's own Cartan dimensions, whose source is the same `child_info`/path algebra), or anything about J. A refutation would be a dropped relation with unchanged dimension (a redundant relation); none occurred because `crels` is minimal here.

## Evidence

Drop-one table (rows per step in `data/`):

| class | steps | true child dimcheck | dropped variants | extra dim +1 / +2 / +3 / +4 | dimcheck rejects | symcheck2 alone says iso |
|---|---|---|---|---|---|---|
| 1 | 16 | 16/16 pass | 86 | 45 / 33 / 7 / 1 | 86 | 86 |
| 2 | 9 | 9/9 pass | 39 | 33 / 6 / 0 / 0 | 39 | 39 |

Step 12 (c1, v = 4, parallel pairs [2, 2], 7 relations): 7 of 7 drops are +2 (one +3), all rejected. End-to-end `compare2` via `toolsmith_endt2_run.py fail`: 25 of 25 `iso`, each in <= 1 s. (Class-2 steps have no parallel arrows, so the c2 control is the plain-quiver case.) The pickles were rebuilt (c1: 16 fail steps, c2: 9, as recorded in E-151/E-174).

## Reproduction

Pickles: `DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 CLS 20000 /tmp/tsm/cCLS.pkl 100`, rerun until it prints "fail steps" (c1 ~ 810 s over two runs, c2 ~ 640 s).

```
CLS=1 timeout 10m .venv/bin/python workshop/rounds/059/toolsmith_dimcheck.py   # ~1 min, asserts, then drop-one table
CLS=2 timeout 10m .venv/bin/python workshop/rounds/059/toolsmith_dimcheck.py   # ~20 s
CLS=1 timeout 10m .venv/bin/python workshop/rounds/057/toolsmith_endt2_run.py fail   # 16/16 iso; CLS=2: 9/9 iso
```

## Prior record

E-177 (dropping a relation accepted 47/47, request for this assertion and for the step-12 / class-2 control); E-174 (25/25 iso). New: the assertion exists and the control is 125/125 including step 12 and class 2. Not in `RETRACTIONS.md`.

## Code changed

`workshop/rounds/057/toolsmith_endt2.py`: added `quotient_dims` (copied from `rounds/058/skeptic_wrongalg.py`, a dead loop removed) and `dimcheck`, called in `compare2` (new key `res['dimcheck']`, verdict string `NO (dim K Q/I' != dim End(T))`). Behaviour on E-174 inputs is unchanged (25/25 iso above). No library (`quivermutation/`) file touched, so no pytest run applies.

## Next

Reasonable now: record the assertion as the soundness condition of E-174/E-177. The n = 10 E-169 edges and the 370 path edges were not rerun with `dimcheck` (cheap: `toolsmith_endt2_run.py paths 0:N`); a skeptic could ask for it. Equal-dim wrong algebras at step 12 (two pairs, a 2-parameter family) are still untested.
