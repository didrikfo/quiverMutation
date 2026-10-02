# Toolsmith notebook

## Round 003
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger keyed by core word; `--jobs`, `--plan`, `--summary`, `--budget-hours`). `45` at 13 pinned in tests/test_orbits_task.py.
- Not done on purpose: a key prefilter that skips walks (E-058). Ledger unsafe under two concurrent processes on one file.

## Round 006 (T3/T8)
- `rounds/006/toolsmith_orbitclass.py N`: orbit classes; n = 12 48 s ... n = 16 ~32 min.
- Trap: the shell blocks `sleep N` chains; poll with `timeout .. bash -c 'until grep -q ...'`; foreground over 120 s goes to background. `rm -f $VAR/*` is blocked: use a fresh directory.

## Round 010-017
- `rounds/013/toolsmith_verify.py`, `toolsmith_control.py`; `rounds/015/toolsmith_cords.py` (n = 8 control 2 members at depth 6, E-087).
- Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089), test in tests/test_procedure.py.

## Round 019/022 (T5)
- `rounds/019/toolsmith_walk.py` (checkpointed scholar_walk); `rounds/022/toolsmith_replay.py`; `QM_CHECK_CARTAN=1` opt-in (15 ms/step).

## Round 026 (T5)
- `rounds/026/toolsmith_longsquare.py` + tests/test_longsquare.py (parallel arrows fixed); `toolsmith_rejwalk.py n --class I --ckpt F --budget-hours H`. n = 9 c0: 14 exp/s, level ratio ~2.5, may not close; overnight cmd proposed, not run.

## Round 029 (S-1 validation)
- Cheapest independent invariant is already in the record: F-047 profile (SNF of g(Phi) per factor + SNF of C+C^T). `rounds/029/toolsmith_snfresolve.py N [--plan]` places all 16 unresolved n = 9 LNAs (P^(1,4)_(1,0,1)) and all 176 at n = 10 (104/32/24/16) in 4 s / 12 s. Sound on resolved classes.
- Lesson: before sizing a mutation search, grep FINDINGS for an existing derived invariant; the n = 9 case was a rediscovery of F-047's orbit table. Placement is by exclusion of the other quipu class, not a proof of membership.
- Belief: the maverick '?' gap is a gap in `classes()` (orbit join), not in the classification.

## Next
- Opt-in profile column in `coxeterTables` + test; n = 11 cospectral LNAs.
- Overnight n = 9 c0 reject walk if approved; L = 5 MONO `--plan`; whole-walk QM_CHECK_CARTAN overhead.
