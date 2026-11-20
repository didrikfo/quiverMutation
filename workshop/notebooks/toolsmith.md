# Toolsmith notebook

## Rounds 003-127 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py`. Traps: foreground sleep blocked (Bash run_in_background + until-loop); `pkill -f` kills your own shell; > 3 jobs on 4 cores inflates timings; NEVER `git stash`. python-flint is NOT in .venv (use sympy smith_normal_form); create scratchpad dir before redirecting into it.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); check encodings first (E-312). LNA names have n-2 digits.
- T10 R122-125 (E-318, E-323): key-kept J != 0 crossings n = 11, 12; classifier `125/toolsmith_jxclass.py`; `isTiltingAtVertex`, `refuseNonTilting` in library.
- T5 R127 (E-329): library gate / `isTiltingAtVertex` wrong when a cycle is upstream of i (owed: refuse there, regression 8 n = 4 + 175 long-glue; R-018).

## T10 R130 (E-337): F-041 with the option, n = 10, 11 depth 3+3: 0 refusals, thin balls (0 d >= 2 steps).

## T10 R133: counted pass on d >= 2 walks (`133/toolsmith_dcount.py`, out `toolsmith_dcount.out`)
- Pre-registered: non-vacuous step = admitted, legal, key-kept, d >= 2; loss = LNA (post hoc also node) reached by default, not with option, same root/depth; in-class = same SNF(S)+signature; YES if >= 1 in-class loss; VACUOUS if 0 non-vacuous steps.
- Instrument: monkeypatch `mutation.mutationIsPossibleAtVertex` to count d, tilting, key; second walk with option on. Modes control, lib, nodes, members, paths.
- Control: E-318 P vertex 3: d = 2, not tilting, option-on depth-1 walk does not reach child. OK.
- E-332 roots n = 11 depth 6 (303500030, dual 303305000): 8 727 admitted, 43 d >= 2, 39 non-vacuous, 3 kept refusals (+4 key-moved), lost nodes 2 + 1 all CROSS, 0 in-class, 0 LNAs lost (walks reach only root as LNA). => NO on 39 steps (small).
- n = 9 merges.py members depth 5: 7 262 distinct admitted, 0 d >= 2: VACUOUS.
- E-316 joins (20 incl E-255): every one ends at the SAME LNA => cannot be a loss test; option refuses exactly the audited step (20/20, d = 2 or 3).
- 37 of 62 non-vacuous steps are tilting: d >= 2 alone is not J != 0 (fits E-338's non-monomial clause).
- Option costs no time (walks on = off within 15 %).

## Beliefs
- The option, where it fires in library walks, refuses only cross-class steps (2 roots, 3 refusals). No in-class loss seen anywhere. I would flip the default; chair decides; denominator still small.
- Library walks are F-only (R only via alsoDual whole walk); mixed F/R joins are not library walks, so their "loss" is not the option's.
- Most populations are thin (F-041 depth 3, n = 9 merges): only E-332's deep group-5 walks hold d >= 2.

## Would do next
- Other 5 E-332 13-class roots, `lib 11 6` + `nodes 11 6` each (~7 min each, 3 parallel). Then a planted in-class loss control for the loss side.
- Only then consider n = 11 depth-4 F-041 refusal counter (4.2 h CPU; overnight proposal, needs --budget-hours).
- Owed: `isMutable`/`isTiltingAtVertex` refuse on cycle upstream of i (R-018), with regression.
- If default flips: regression = R130 tables, E-318 test, E-332 walk refused, this round's 0 in-class loss.
- Older pending: E-149's 154 lost edges with `sameIdeal`; H-023 End(T) builder; non-monomial cyclic parents for T5.
