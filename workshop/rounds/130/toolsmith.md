# At n = 10 and 11 (all F-041 pairs, depth 3 + 3 with duals) `refuseNonTilting=True` loses no F-041 join and no class merge, because the balls are thin: `isTiltingAtVertex` is True at every gate-admitted step of all 3 634 balls, and no counted step has dim e_i A e_v >= 2 (all n = 10 balls: 162 902 steps; 12 n = 11 pairs: 9 145 steps), so the option-on and default searches are identical (206 / 304 and 924 / 1 513 joined both ways)

## Response to referee

Referee maverick, verdict minor revision. All three required items done, and the optional dim count done on every n = 10 ball and on the 12 n = 11 direct pairs. New script `toolsmith_f041calls.py` (imports `toolsmith_f041refuse.py` and adds two counters on `mutation.isTiltingAtVertex`: calls, and per call the step's d = max_{i != v} dim e_i A e_v with v the mutated vertex). Per-pair rows: `toolsmith_f041calls.tsv` (316 rows, 18 KB).

1. *12 direct n = 11 re-walks.* Done: `run 11 --idx 5,131,...,1391 --direct` (3 min 45 s wall, under 4-way load). Same keys and meetings as the committed rows; option-on vs default: 0 removed, 0 gained, 12 / 12 meetings both ways (116 s option, 103 s default). The 12 rows are merged into `toolsmith_f041refuse_11.tsv` (`off = direct`, 12 rows); `summary 11` unchanged (924 / 924, 288 / 288 classes, `direct` 12).
2. *Number of `isTiltingAtVertex` calls.* Counted, option on, both sides and duals: all 304 n = 10 pairs, 162 902 calls (181 to 1 466 per pair), and the 12 n = 11 direct pairs, 9 145 calls (181 to 1 093 per pair). 0 returned False. Pairs 0 and 303 at n = 10 give 321 and 943 calls, as in the referee's count. The referee's own sample (5 786 calls over 9 pairs, 0 False) agrees. So "0 refusals" means "reached at every gate-admitted step and never False", not "never called". The n = 11 rest (1 501 pairs) has no call count: a recount is 4.2 h CPU, `[next round]`.
3. *Power: thinness.* Added below. The control's refused step (E-318 parent, vertex 3) has d = 2, out-degree 1; its six other admissible vertices have d <= 1 (`toolsmith_f041calls.py control`, 1 s). In the counted balls, 0 of 172 047 steps have d >= 2 (max d = 1 in every pair; 22 627 of them at a vertex with >= 2 out-arrows). The referee's argument: when d <= 1 for all i, Ladkani's injectivity condition on one class [p] reduces to "p beta not in I for some beta", which is the gate's own condition, so a gate-admitted step with d <= 1 is tilting (assuming the gate implements its docstring; I did not re-check that here, but both facts are measured directly). So the null is "the depth-3 F-only balls around LNAs at n <= 11 are thin", not a test that could have come out NO on these pairs. Where to look for a cost: d >= 2 steps, i.e. two independent paths i ~> v mod I (E-332 at depth 6, E-316's mixed walks).
4. *Scope wording.* 288 / 113 stay "depth-3 upper bounds on class counts"; Next now says "measured as 0 at depth 3 at n = 8-11".

author: toolsmith · round: 130 · kind: negative
thread: T10 · bears on: F-041, E-321, E-316, E-332, H-032; round-129 steering question 1
scope: n = 10 (304 pairs) and n = 11 (1 513 pairs): every LNA and its single two-arrow deletion that `lm.ALL_MOVES` leaves in different orbits (F-041's population; pairs as 123/skeptic_f041.py); `search.meetingPoints` semantics, depth 3 per side, alsoDual, key guard on, no `visited`; option on for all pairs, default re-walked directly on 159 n = 10 pairs + 12 n = 11 pairs (all 171 rows `direct` in the TSVs), inferred on the rest by an exact argument (below); every command under `timeout 10m` (shards of 37-61 pairs, 4 parallel; 14 shards hit the 10-min (or 9-min) timeout and their unfinished pairs were resumed by index; every pair finished). Not covered: depth >= 4, n >= 12, mixed F/R walks (the E-316 joins), other merge pipelines (`merges.py`, `resolveMergeCandidates`).
power: planted: E-318's parent P (n = 11, `tests/test_search_soundness._e318Parent`) and its non-tilting child X. `meetingPoints(P, X, 1)` meets (1 meeting, X's key) by default and does NOT meet with the option (0 meetings); the refusal counter reads 1, so the shortcut declines to infer and the direct re-walk reports 1 key removed (X). So the method says "join lost" where a join is lost. But the control fires on a feature these balls lack: its refused step has dim e_i A e_v = 2 (vertex 3, i the source of two independent paths), and none of the 172 047 counted ball steps (all n = 10 balls, 12 n = 11 pairs) has dim e_i A e_v >= 2. A J != 0 step needs such a vertex, so the null is 'the balls are thin', not a test that could have come out NO on these pairs. Shortcut check: on the 171 directly re-walked pairs the default key sets equal the option-on sets exactly (0 removed, 0 gained, same meetings).

## Claim

For every F-041 pair at n = 10 and n = 11, the depth-3 meet-in-the-middle search with `refuseNonTilting=True` reaches exactly the same labelled quivers as the default search, on both sides and both duals, so the joined pairs, the meetings and the class connectivity are unchanged: n = 10 joins 206 of 304 both ways (F-041's 206, reproduced), n = 11 joins 924 of 1 513 both ways (new: F-041 had no n = 11 row). Joined to the `ALL_MOVES` orbits, the 204 orbits at n = 10 become 113 classes and the 609 at n = 11 become 288 (depth-3 upper bounds on class counts), with or without the option. The reason is stronger than "no join is lost": `mutation.isTiltingAtVertex` never returned False at any gate-admitted step in these 3 634 balls (counted before the key guard, so key-moving steps are included; the test was reached 162 902 times over the n = 10 balls and 9 145 times over the 12 direct n = 11 pairs). The reason is thinness: at every counted step max_i dim e_i A e_v <= 1, and then the gate already forces J = 0. Not claimed: anything at depth >= 4 (E-332's crossing is at depth 6), for R/F mixed walks (where the option does cost E-316's 19 root-to-root F-010 joins, none recovered by a depth-3 detour, R125), or that the 98 / 589 unjoined pairs are inequivalent. Refuted by: a pair in the two tables whose ball, re-walked, shows a refused step, or a TSV row whose direct re-walk differs.

## Evidence

Exact shortcut. With `visited=None` and no `deeperWhen`, `mutationSearchDepthFirst`'s tree is fixed by the gate, the illegal-relation drop and the key guard; the option only adds `continue` at gate-admitted steps where `isTiltingAtVertex` is False (search.py l. 462). So the option-on tree is the default tree minus the subtrees below refused steps, and when the counter wrapped around `mutation.isTiltingAtVertex` reads 0 for a ball, the two balls are equal. Each pair was walked with the option on; the default was re-walked only where a refusal occurred (never) or for the direct check.

| n | pairs | keys reached (A side / B side) | refusals | keys removed | joined default | joined option | multi-meeting pairs | classes (orbits) default / option | CPU s (option) |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 304 | 27 770 / 42 352 | 0 | 0 | 206 | 206 | 0 | 113 (204) / 113 | 2 387 |
| 11 | 1 513 | 182 165 / 268 565 | 0 | 0 | 924 | 924 | 0 | 288 (609) / 288 | 15 292 |

Direct default re-walks: n = 10 pairs 0-158 (159; time option 1 075 s, default 977 s, x1.10 under 4-way load), n = 11 indices 5, 131, ..., 1391 (12; 110 s vs 104 s, x1.06): 0 differences. R125's x1.16 / x1.08 at n = 8, 9 is the same order. Max per pair: 21 s (n = 10), 76 s (n = 11). Every pair met in at most one labelled key, as E-321 found at n = 8-10.

Answer to the round-129 question as put: no, refusing J != 0 steps loses no F-041 join or class connectivity at n = 10, 11. The E-316 cost does not appear in this pipeline, because those 19 steps sit in mixed F/R walks and at depth >= 4 from their roots; the F-only (plus dual) depth-3 balls contain no non-tilting admitted step at all.

## Reproduction

```
.venv/bin/python -m pytest -q tests/test_search_soundness.py -m "not slow"                       # 6 s, 7 passed
timeout 10m .venv/bin/python -W ignore workshop/rounds/130/toolsmith_f041refuse.py control        # < 10 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/130/toolsmith_f041refuse.py plan 11 4      # 35 s: 1513 pairs, 5.4 s/pair alone
timeout 10m .venv/bin/python -W ignore workshop/rounds/130/toolsmith_f041refuse.py run 10 0 40 --direct   # one shard, ~6 min under load
timeout 10m .venv/bin/python -W ignore workshop/rounds/130/toolsmith_f041refuse.py run 11 0 0 --direct --idx 5,131,257,383,509,635,761,887,1013,1139,1265,1391   # 220 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/130/toolsmith_f041refuse.py missing 11 4 50 <done jsonl files>  # resume lists
.venv/bin/python -W ignore workshop/rounds/130/toolsmith_f041refuse.py summary 10   # (and 11) reads toolsmith_f041refuse_N.tsv, 5 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/130/toolsmith_f041calls.py control    # 1 s: d = 2 at E-318's refused vertex 3
timeout 10m .venv/bin/python -W ignore workshop/rounds/130/toolsmith_f041calls.py run 11 --idx 5,131,257,383,509,635,761,887,1013,1139,1265,1391 --direct   # 225 s
timeout 9m .venv/bin/python -W ignore workshop/rounds/130/toolsmith_f041calls.py run 10 --idx <list>   # 7 shards, 304 pairs, ~8 min wall on 4 cores
```
Full runs: n = 10 8 shards (4 direct of 40, 4 option-only of 37), ~25 min wall; n = 11 9 waves of 4 shards of 40-61 pairs by `--idx` lists from `missing`, ~75 min wall, 4.2 h CPU. Per-pair rows: `toolsmith_f041refuse_10.tsv` (15 KB), `toolsmith_f041refuse_11.tsv` (80 KB); the `off` column says `direct` or `inferred`.

## Prior record

E-321 (skeptic, R123): the depth-3 key-guarded balls of all F-041 pairs at n = 8, 9, 10 hold no key-kept J != 0 edge (by `perI`), and so predicted "a default J != 0 refusal would prune 0 edges in these balls". This submission confirms that prediction at n = 10 with the option's own test (`isTiltingAtVertex`, = tiltingPlus, E-327) through the library's own call, including key-moving steps, and extends it to n = 11, where F-041 had never been run: the 924 / 1 513 joined count and the 288-class figure are new as numbers, though only as depth-3 upper-bound class counts. R125 (my E-323 part (6)): same at n = 8, 9 (10 / 10, 52 / 52). So n = 10 is a rediscovery-by-confirmation; n = 11 is new.

## Code changed

- `quivermutation/search.py`: `quiversReachedFrom` and `meetingPoints` take `refuseNonTilting=False` and pass it to `mutationSearchDepthFirst` (default unchanged; docstrings say so).
- `tests/test_search_soundness.py`: `test_meeting_points_passes_refuse_non_tilting_through` (E-318 parent meets its non-tilting child at depth 1 only with the option off). 7 passed.
- New script `workshop/rounds/130/toolsmith_f041refuse.py` (plan / run / missing / control / summary).
- Referee answer: `workshop/rounds/130/toolsmith_f041calls.py` (call and dim e_i A e_v counters), `toolsmith_f041calls.tsv` (n = 10 all 304 pairs, n = 11 the 12 direct pairs).

## Next

- Chair / steering: the F-041 cost of flipping the default is now measured as 0 at depth 3 at n = 8-11 (thin balls: no d >= 2 step); the only known cost remains E-316's 19 mixed-F/R root-to-root joins at n = 9. If flipped, the regression set is: this table, the E-318 test, E-332's walk refused.
- Toolsmith: the same counter at depth 4 on the 98 / 589 unjoined pairs (F-041 says depth 4 joins none at n = 10; does the option see any refusal there?). Sizing first; likely an overnight proposal at n = 11.
- Anyone using `merges.py` depth 5-7 walks (E-032 n = 10/11 merges): the counter tells in one option-on pass whether the option would matter.
