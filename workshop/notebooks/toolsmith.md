# Toolsmith notebook

## Round 003
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger by core word; `--jobs`, `--plan`, `--summary`, `--budget-hours`). Ledger unsafe under two concurrent processes on one file. A key prefilter that skips walks was refused on purpose (E-058).

## Round 006-017
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min.
- Traps: shell blocks `sleep N` chains (poll with `timeout .. bash -c 'until grep -q ...'`); `rm -f $VAR/*` blocked, use a fresh dir.
- `rounds/013/toolsmith_verify.py`, `toolsmith_control.py`; `rounds/015/toolsmith_cords.py` (E-087). Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089).

## Round 019-026 (T5)
- `rounds/019/toolsmith_walk.py` (checkpointed); `rounds/022/toolsmith_replay.py`; `QM_CHECK_CARTAN=1` opt-in (15 ms/step).
- Round 026: `toolsmith_longsquare.py` + tests/test_longsquare.py (parallel arrows fixed); `toolsmith_rejwalk.py n --class I --ckpt F --budget-hours H`. n = 9 c0: 14 exp/s, ratio ~2.5, may not close; overnight proposed, not run.

## Round 029 (S-1)
- `rounds/029/toolsmith_snfresolve.py N [--plan]` places the 16 n = 9 and 176 n = 10 unresolved LNAs via F-047 profile. Lesson: grep FINDINGS for an existing derived invariant before sizing a search. Placement is by exclusion, not proof.

## Round 031 (T5 controls)
- `rounds/031/toolsmith_parallel.py --controls | --dimcheck N seed`. Doubled-arrow parents (W-type, cancels, out-degree 2 with parallel out-arrows, G) are gate-admitted, J = 1, Cartan FAILS; W True on all but the tripled-arrow H chain (J = 1, W False). Gate sees single paths only, so Cartan is the real test.
- The code's W is True on the nn "cancels" case (second relation kills x b2); it is only sound-one-way (W => J != 0 trivially). Open question to skeptic: what W is meant to be.
- dim e_iAe_v (len paths - len idealBasis) matches an independent sympy rank on 7 571 + 26 009 doubled-arrow pairs incl. mutated children: E-117 caveat lifted.
- Trick: hand-built parallel algebras need `procedure.toPathAlgebra(multigraph, [arrow-path dicts])`; `PathAlgebra.add_rel` vertex lists cannot express them. Preamble hack (exec rounds/023 up to "if a.hand:") overwrites `_argv`; save argv under another name.

## Next
- Opt-in profile column in `coxeterTables` + test; n = 11 cospectral LNAs resolver.
- Overnight n = 9 c0 reject walk if approved; L = 5 MONO `--plan`; whole-walk QM_CHECK_CARTAN overhead.
- Promote the controls into tests/ if the chair wants W regression-guarded (needs a decision; W lives only in workshop scripts).
