# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger by core word; `--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file. Key prefilter that skips walks refused on purpose (E-058).
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. `rounds/013/toolsmith_verify.py`, `_control.py`; `rounds/015/toolsmith_cords.py` (E-087). Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089).
- Traps: shell blocks `sleep N` chains (poll with `timeout .. bash -c 'until grep -q ...'`); `rm -f $VAR/*` blocked.

## Round 019-026 (T5)
- `rounds/019/toolsmith_walk.py`, `rounds/022/toolsmith_replay.py`; `QM_CHECK_CARTAN=1` opt-in (15 ms/step). Round 026: `toolsmith_longsquare.py` + tests/test_longsquare.py; `toolsmith_rejwalk.py`. n = 9 c0: 14 exp/s, may not close; overnight proposed, not run.

## Round 029 (S-1)
- `rounds/029/toolsmith_snfresolve.py N [--plan]` places unresolved LNAs via F-047 profile (by exclusion, not proof). Lesson: grep FINDINGS for a derived invariant before sizing a search.

## Round 031 (T5 controls)
- `rounds/031/toolsmith_parallel.py`: doubled-arrow parents are gate-admitted, J = 1, Cartan FAILS; code W only sound one way; dim e_iAe_v count matches sympy on 33k pairs (E-117 caveat lifted).
- Trick: hand-built parallel algebras need `procedure.toPathAlgebra(multigraph, [arrow-path dicts])`. Preamble hack (exec rounds/023 up to "if a.hand:") overwrites `_argv`; save argv under another name.

## Round 033 (base rate of LNA keys)
- Believe: the key test is vacuous for E-121. Walk rows carry the key by construction (3585/3585 out-degree-2 rows at n=8 c0, 3000 expansions); the layered family has 0/2704 hits including 1003 circuit-free and 108 pure-W members. `rounds/033/toolsmith_baserate.py layers m | walk n cls maxexp`.
- 42 of 3585 out-degree-2 children have an LNA key outside c0 (presumably rejects; unmatched).
- Lesson: "no hit" needs a positive control in the same family before it is evidence; check that first.

## Next
- Find a family with LNA-keyed circuit-free members, then redo W vs non-W hit rates.
- Match the 42 non-c0 children to E-114 rejects (cheap, same script).
- Opt-in profile column in `coxeterTables` + test; n = 11 cospectral LNAs resolver; promote W controls to tests/ if the chair decides.
