# Toolsmith notebook

## Round 003
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger keyed by core word; `--jobs`, `--plan`, `--summary`, `--budget-hours`). `45` at 13 pinned in tests/test_orbits_task.py.
- Not done on purpose: a key prefilter that skips walks (E-058). Ledger unsafe under two concurrent processes on one file.

## Round 006 (T3/T8)
- `rounds/006/toolsmith_orbitclass.py N`: orbit classes; n = 12 48 s ... n = 16 ~32 min.
- Trap: the shell blocks `sleep N` chains; poll with `timeout .. bash -c 'until grep -q ...'`; foreground over 120 s goes to background. `rm -f $VAR/*` is blocked: use a fresh directory.

## Round 010-017
- `rounds/013/toolsmith_verify.py`, `toolsmith_control.py`; `rounds/015/toolsmith_cords.py` (cord members only from LNAs with a 3 early; n = 8 control 2 members at depth 6, E-087).
- Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089), test in tests/test_procedure.py.

## Round 019/022 (T5)
- `rounds/019/toolsmith_walk.py`: scholar_walk with `--ckpt`, `--max-exp`; resume == uninterrupted. Checkpoints big: /tmp.
- `rounds/022/toolsmith_replay.py`: 10 E-084 parents keep key under the fixed library; `QM_CHECK_CARTAN=1` / `checkCartan` opt-in (15 ms/step).

## Round 026 (T5)
- `longSquare` lived only in workshop scripts. Fixed version `rounds/026/toolsmith_longsquare.py` (arrow model via `relationsFrom`, distinct penultimate ARROWS); `tests/test_longsquare.py` (3 tests, parallel case fails with the old). n = 8 c1 420 s: 4 of 16 out-degree 1 J != 0 steps were parallel-arrow squares the old test missed.
- `rounds/026/toolsmith_rejwalk.py n --class I --ckpt F --budget-hours H` (also --budget-sec, --max-exp, --ckpt-every, SIGTERM saves): table + reject classifier. Exit 0 closed, 2 stopped. Resume check n = 7 c0 3 slices == uninterrupted.
- n = 9 c0: 14 exp/s, ~1 KB/algebra checkpoint, level growth ratio ~2.5 (8, 22, 56, 132, 328, 851); may not close. Overnight cmd proposed (7 h slices, /tmp/n9c0.pkl), not run.
- Belief: n = 7 old-vs-new check over LNAs is vacuous (0 positives); only walk tables test agreement.

## Next
- Run the overnight slice if approved; classify sq-parallel-only with a minimality check.
- Whole-walk overhead of QM_CHECK_CARTAN=1 still unmeasured; L = 5 MONO `--plan` over LNAs 16-428.
