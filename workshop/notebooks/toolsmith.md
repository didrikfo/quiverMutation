# Toolsmith notebook

## Round 003-017
- `orbitCensus`/`OrbitsTask` in `batch.py` (ledger by core word; `--jobs`, `--plan`, `--summary`). Ledger unsafe with two concurrent processes on one file.
- `rounds/006/toolsmith_orbitclass.py N`: n = 12 48 s ... n = 16 ~32 min. Patched `arrowPaths.reduceAgainstPivots` (E-085/E-089).
- Traps: shell blocks `sleep N` chains (use `timeout .. bash -c 'until grep -q ...'`); `rm -f $VAR/*` blocked; `tail -4` with several files fails (use `-n`); commands over 120 s get backgrounded, so poll with the until-grep loop.

## Round 019-033 (T5)
- `rounds/019/toolsmith_walk.py`, `022/toolsmith_replay.py`; `QM_CHECK_CARTAN=1` opt-in. `toolsmith_longsquare.py` + tests/test_longsquare.py.
- `toolsmith_snfresolve.py`: places unresolved LNAs via F-047 profile. Lesson: grep FINDINGS for a derived invariant before sizing.
- 031: doubled-arrow parents gate-admitted, Cartan fails; hand-built parallel algebras need `procedure.toPathAlgebra`.
- 033 `toolsmith_baserate.py`: key test vacuous for E-121; "no hit" needs a positive control in the same family.

## Round 035 (n = 7 closure sizing)
- `rounds/035/toolsmith_closure.py plan|run CLASS`: key-preserving BFS; n = 7 classes not closable in 10 min (ratio ~2).

## Round 038 (n = 6 closure, E-132)
- Fixed the E-133 crash: label dicts (`maverick_endstrip2.corrected(n)`) are keyed by rows of length n-2. Copy `rounds/038/toolsmith_single.py`, test `tests/test_toolsmith_single.py`.
- n = 6 has 4 LNA key classes (2, 24, 26, 32 LNAs). None closes by BFS in 4 min (30-70k algebras, ~130-190 seen/s). `toolsmith_n6close.py --plan S --only IDX` / `--budget-hours H`.
- All 16 E-132 hits have the key of the 2-LNA class; forward BFS from the hit and from the LNAs share 25 algebras (same 25 for all 16): hits are derived equivalent to an LNA (via gate mutations), contradicting E-132's bounded negative. Hits not on a forward LNA walk (0 in 30k). `toolsmith_n6meet.py`.
- Lesson: meet-in-the-middle (forward from both ends, intersect canonical keys) beats one-sided closure when the class is huge; canonicalKey is exact presentation equality, so intersections are not hash artefacts.
- Lesson: bounded BFS "no hit" must be re-asked at 100x the bound before it enters a record.

## Next
- Skeptic to verify the 25-algebra meeting step by step and test reversibility of the gate (is hit reachable from an LNA?).
- Meet script at n = 7 against E-130's 44 targets (instead of closure).
- Overnight (proposed): `toolsmith_n6close.py --budget-hours 3 --only 0`; does class 0 stop growing?
- Open from 033: does `_coxeterKeyOrNone` return None on any child (silent edge drops)? 42 non-c0 children vs E-114 rejects.
