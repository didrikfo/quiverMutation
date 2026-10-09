# The n = 8 class 2 depth-8 walk completes under the fix: 0 key-moved steps; the 10 of E-086 became guard-admitted tilting steps

author: toolsmith · round: 019 · kind: tool + result
thread: T5 · bears on: E-086, E-087, E-091, E-092, F-052 (as the caution)

## Claim

`workshop/rounds/019/toolsmith_walk.py` is `workshop/rounds/014/scholar_walk.py` with checkpoint/resume (`--ckpt FILE`,
atomic pickle of seen-set, current/next frontier with paths, position in the level, counters; exit 2 when the slice
budget is spent). Sliced walks give output byte-identical to the uninterrupted walk (n = 7 classes 0, 1, depth 6; 6 and 10 slices).
With it, depth 8 of n = 8 class 2 under the library as it now is (E-091 fix) ran to the end of the level in two slices
(493 s + 208 s, one process): 24 316 expansions, 63 221 distinct algebras, 104 629 guard/tilt steps, 2 rejections
(gate True, `tiltingPlus` False), **0 "gate+tilt but key moves"**. Over the same first 20 899 expansions that E-086
completed (E-086's 560 s cap), the counts are: guard-tilt 89 189 now vs 89 179 + 10 key-moved in E-086
(total 89 189 both); noguard-NOTtilt 1 and the rejection line identical; distinct algebras 54 333 vs 54 326.
So the 10 key-moved steps are exactly the steps that now keep the key: they were the E-087 defect, as E-087's replay said,
and this is now seen by the walk itself rather than a replay of 10 saved parents.
Not claimed: the walk is not closed (frontier 38 907 at depth 8; depth 9 not run); the second rejection
(expansion > 20 899) was not replayed or inspected; a plausible-but-untested reading of the "+7" in distinct algebras is
that the 10 fixed children are new nodes or collapse into seen ones (3 of 10 collapsed); I did not check which.

## Evidence

| run | expansions | distinct | guard,tilt | noguard,tilt (key moved) | rejections |
|---|---|---|---|---|---|
| E-086 file (library before E-091), cap 560 s | 20 899 | 54 326 | 89 179 | 10 | 1 |
| this, fixed library, `--max-exp 20899` | 20 899 | 54 333 | 89 189 | 0 | 1 |
| E-092 (monkeypatch), budget 585 s | 17 058 | 43 638 | 71 659 | 0 | 1 |
| this, fixed library, depth 8 complete | 24 316 | 63 221 | 104 629 | 0 | 2 |

Resume check (library as it is, `--depth 6`, `--budget-sec 3`): n = 7 class 0 (6 slices) and class 1 (10 slices):
the summary (counts, rejections with their path tuples) equals `scholar_walk.py`'s output with `diff`
(`toolsmith_walk_resume_check.txt`). The `--max-exp` run is itself a resume (3 slices), and its depth-1..7 counts equal E-086's lines.
The second rejection: parent path (17, 8, 5, 6, 8, 8, 2, 5), vertex 5, in `toolsmith_walk_n8_c2_slice2.txt`.
Checkpoint at depth 8 is 38 MB: kept in /tmp (`/tmp/tw_n8c2.ckpt`), not in the repository.
Timing: depth 8 costs 235 s to depth 7 then 440 s for the level; a slice of 480-500 s fits the 10-minute limit.

## Reproduction

```
# resume equals uninterrupted (about 2 min each class):
timeout 10m .venv/bin/python workshop/rounds/014/scholar_walk.py 7 --class 0 --depth 6 > /tmp/a.txt
rm -f /tmp/c.ckpt; until timeout 10m .venv/bin/python workshop/rounds/019/toolsmith_walk.py 7 --class 0 --depth 6 --budget-sec 3 --ckpt /tmp/c.ckpt > /tmp/b.txt; ! grep -q CHECKPOINT /tmp/b.txt; do :; done
diff <(grep -v '^depth' /tmp/a.txt) <(grep -v '^depth\|RESUME' /tmp/b.txt)
# n = 8 class 2, depth 8, two slices (493 s, 208 s):
timeout 10m .venv/bin/python workshop/rounds/019/toolsmith_walk.py 8 --class 2 --depth 8 --budget-sec 480 --ckpt /tmp/tw_n8c2.ckpt   # run twice
# same prefix as E-086 (3 slices of <= 420 s): add --max-exp 20899 and repeat until the summary prints
```
Outputs: `toolsmith_walk_n8_c2_slice1.txt`, `..._slice2.txt`, `..._first20899.txt`, `toolsmith_walk_resume_check.txt`.

## Prior record

E-086 (the 10 key-moved steps, truncated depth 8), E-087 (defect of `reduceAgainstPivots`, replay of 10 parents),
E-091 (fix), E-092 (re-run stopped at 17 058, "neither reproduced nor refuted"). This closes E-092's gap: refuted under the fix,
for the whole first 20 899 expansions and the remaining 3 417 of depth 8. Not new in kind: it confirms E-087's diagnosis.
Correction to the record's wording: E-086's n = 8 class 2 depth 8 was itself partial (20 899 of 24 316 expansions), so
"20 899 depth-8 expansions" in E-092 is cumulative expansions in a truncated level, not the level size.

## Code changed

New file only: `workshop/rounds/019/toolsmith_walk.py`. Library untouched, no tests run (none touch it; check is the diff above).
Behaviour differences from `scholar_walk.py`: algebra paths held in frontier tuples instead of `id()`; the summary prints also on a checkpointed slice;
budget is per slice. Known limit: a checkpoint is a pickle of library objects, valid only for the same code revision.

## Next

Experimentalist/chair: n = 8 class 2 depth 9 (frontier 38 907; est. 3x depth 8, about 35 min, 4+ slices) as an `OVERNIGHT.md` item
with `--budget-hours` (the script has only per-slice `--budget-sec`; a shell loop of slices is the proposal). Scholar: inspect the
second rejection (parent path (17, 8, 5, 6, 8, 8, 2, 5)) against E-086's A5 shape. Toolsmith: add `--budget-hours` wrapper if wanted.
