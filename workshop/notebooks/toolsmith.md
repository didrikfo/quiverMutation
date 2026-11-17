# Toolsmith notebook

## Rounds 003-125 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py`. Traps: foreground sleep blocked (Bash run_in_background, then wait for the notice); `pkill -f` kills your own shell; > 3 jobs on 4 cores inflates timings (~1.6x); check `ps` after a timeout. NEVER `git stash`.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); check encodings before tracing code (E-312). LNA names have n-2 digits.
- T10 (R122, R125; E-318, E-323): key-kept J != 0 crossings at n = 11 (13-class) and n = 12 (47-group). Classifier `125/toolsmith_jxclass.py`. Library R125: `isTiltingAtVertex`; `refuseNonTilting` opt-in in `mutationSearchDepthFirst` / `linesReachedFrom`.

## T5 / H-023, round 127 (`127/toolsmith_cycgate.py`)
- Loops at i never admitted; walks never step from cyclic parents. 12 136 admitted cyclic-monomial steps checked; D(i)-letter lemma. Library gate and `isTiltingAtVertex` wrong when a cycle is upstream of i (8 n = 4 cases, 175 long-glue); warning added; neither reaches a walk. Monomial: exact gate == J = 0.

## T10, round 130: F-041 with refuseNonTilting (`130/toolsmith_f041refuse.py`)
- Library: `quiversReachedFrom` / `meetingPoints` now take `refuseNonTilting` (test added in test_search_soundness).
- Exact shortcut: no `visited`, no `deeperWhen` => option-on tree = default tree minus subtrees under refused steps. So wrap `mutation.isTiltingAtVertex` with a counter, walk ONCE with the option on; counter 0 => default identical. Halves cost. Verified on 171 direct re-walks, and the E-318 planted meeting (counter 1, join lost).
- Result: n = 10 (304 pairs) and n = 11 (1 513 pairs), depth 3+3 with duals: 0 refusals in all 3 634 balls (even before the key guard). Joined 206 / 924 both ways; classes 204 -> 113 and 609 -> 288 both ways. n = 10 confirms E-321's prediction; n = 11 F-041 row is new.
- Referee (maverick, minor): merged 12 direct n = 11 rows; counted calls (n = 10 all: 162 902; n = 11 12 pairs: 9 145; 0 False) and d = max_i dim e_i A e_v per step: 0 steps with d >= 2 (E-318 control: d = 2). Thinness: d <= 1 + gate => J = 0. So null = 'balls are thin'. `130/toolsmith_f041calls.py`.
- Lesson: a 0-refusal count needs a call count AND a reason it could have fired (feature present in the population).
- Cost: n = 11 4.2 h CPU, ~75 min wall in 9 waves; pairs vary 4-76 s. Resume by `missing` -> `--idx` lists (interleaved shards balance load). Outputs as TSV (jsonl at n = 11 was 294 KB > 200 KB limit).

## Beliefs
- refuseNonTilting costs nothing measurable in F-only depth-3 balls at n = 8-11 (they are thin: no d >= 2 step); its only known cost is E-316's 19 mixed-F/R root-to-root F-010 joins (n = 9). The crossings it removes are deep (depth 4 with R, depth 6 pure F, E-332). I would flip the default; the human/chair decides.
- Pure-F library walks cross at depth 6 (E-332); depth-3 balls never see a non-tilting admitted step at n <= 11.
- Cyclic/loop part of T5 reduces to acyclic D(i). Remaining T5 gaps: analytic (H3), Hom/End(T), non-monomial parents.
- Split disagreement counts by admitted vs refused and by whether caps allow a YES.

## Would do next
- Look for d >= 2 steps first (cheap): E-332 depth 6, E-316 mixed walks; n = 11 call/d count on the other 1 501 pairs (4.2 h CPU).
- Counter at depth 4 on the 98 / 589 unjoined F-041 pairs (size first; n = 11 likely overnight, script needs --budget-hours, exit 2).
- Counter over `merges.py` depth 5-7 walks (E-032 n = 10/11 merges): one option-on pass says whether the option matters there.
- If the default flips: regression = R130 tables, E-318 test, E-332 walk refused.
- Pending from T10: overnight proposal (55 remaining 59-class roots depth 4; 47-group's 14-class; other n = 12 groups); depth-7 n = 11 13-class; E-149's 154 lost edges with `sameIdeal`; H-023 End(T)-presentation builder.
- Non-monomial cyclic parents for T5; exact gate on cycles (regression: 8 n = 4 + 175 long-glue cases).
