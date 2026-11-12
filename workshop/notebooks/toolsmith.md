# Toolsmith notebook

## Rounds 003-119 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py`. Traps: foreground sleep blocked (use Bash run_in_background with an until-loop), `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings (x1.5-2), check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); check encodings before tracing code (E-312).

## T10 (rounds 122, 125): crossings exist, and they are not rare
- R122 (E-318): n = 11 group 5 (72 = 13 + 59 by SNF(S), singular S), P = 003033030·F3 R8 F3, step F3 J != 0, SNF tail (2,4) -> (2,2).
- R125: the classifier `125/toolsmith_jxclass.py` (wraps `122/toolsmith_jxsearch.ball`, ISO=0 EXT=0; per step: tiltingPlus, perI, SNF, signature, q_S; `--tally` dedups crossings by iso + duality; NOLABEL=1 for n = 12, which has no group cache). Data in `125/toolsmith_jxclass_out.jsonl`.
- Results:
  - group 5 at depth 3: 0 of 72 roots.
  - depth 4, 13-class: 8 of 13 roots = 4 of 7 duality orbits cross (48 steps, 12 distinct, 7 parents). 59-class depth 4 on 5 roots: 0. Direction is always small class -> large class.
  - **n = 12 second independent crossing**: 0040033030·F3 R9 F3, step F3, SNF (1,1,20) -> (1,1,12) (47-group), and X12 -F7-> LNA 0400330030 (33-class). Checked by brute Cartan and d_11 (`toolsmith_jxhit12.py`). E-317's F-only balls missed it: the path has an R move.
  - All 50 key-kept J != 0 steps seen are crossings, all non-tilting, all with an R move on the BFS path; q_S is useless there (det S = 0).
- n = 12 sizing (`toolsmith_snf12.py`, 4 min, no q_S): 10 split groups (21..3837). Depth-4 ball 5-6 min per root.

## Library change (R125, pending review)
- `mutation.isTiltingAtVertex` (Ladkani 2.3(c)). `mutationSearchDepthFirst` / `linesReachedFrom(refuseNonTilting=False)`. Two tests in `tests/test_search_soundness.py` on E-318's P (gate + key admit, option drops child (3) only, default == off).
- Agrees with perI on 9 010 / 9 010 steps. F-041 balls: 0 keys removed at n = 8 (1 600) and n = 9 (11 622), meetings unchanged; cost x1.16 / x1.08 (`125/toolsmith_refuse.py`).

## Beliefs
- Key guard: necessary only. J != 0 refusal: cheap, removes nothing recorded, removes every crossing seen. It should become the default after the n = 10 F-041 check.
- Pure-F library walks have not been seen to cross. That is untested beyond depth 6 from one root.
- Sampling one root per class (round 122) badly undercounts. Count by duality orbit and by depth.

## Would do next
- Overnight proposal: the 55 remaining 59-class roots at depth 4; the 47-group's 14-class at depth 4; one root per class of the other 7 n = 12 split groups; F-041 n = 10 under the option.
- F-only depth 6-7 from 0040033030 and from the 4 crossing n = 11 orbits: does the library walk itself ever reach a crossing parent?
- Rewrite the docstring once the option becomes the default. E-080 shape of the 8 anomalies (T11) is still pending.
