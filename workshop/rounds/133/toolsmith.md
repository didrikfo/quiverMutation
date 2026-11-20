# refuseNonTilting loses no in-class join or node on the one population with refusals (E-332's two n = 11 depth-6 library walks: 39 non-vacuous steps, 3 refusals, all 3 removed nodes cross-class); the n = 9 merges walks are vacuous (0 of 7 262 steps with d >= 2)

author: toolsmith · round: 133 · kind: negative
thread: T10 · bears on: H-032, H-015, F-057, E-332, E-316, E-337, E-338
scope: n = 11 library F-only walks (defaults, `fingerprint.Visited`) from 303500030 and dual(303305000), depth 6; n = 9 `merges.py` leftover members (9 members x member/relation dual = 9 distinct walks), depth 5; recorded mixed F/R paths: 20 E-316/E-255 joins (E-321 replay lists) + E-318 5-step + E-332 10-step walk, path level only. Not covered: other E-332 roots (5 of 7), depth >= 7, mixed F/R BFS balls with the option (the library has no mixed walk), n = 10-12 merges, the n = 11 depth-4 F-041 run (not started, as instructed).
power: refusal side: E-318's planted step (`_e318Parent`, vertex 3) has d = 2, `isTiltingAtVertex` False, and the option-on depth-1 library walk from P does not reach its child (default does): caught. Loss-classifier side: it labels 3 of 3 known E-332 crossing nodes CROSS. No planted IN-class loss was built, so a NO on "in-class loss" rests on the 3 refusals observed, not on a control.

## Claim

On the walks that actually hold steps where the option can fire (d = max_i dim e_i A e_v >= 2, admitted, legal, key-kept: "non-vacuous"), `refuseNonTilting=True` refuses only steps that leave the SNF(C + C^T) class: at E-332's two crossing roots (n = 11, depth 6) there are 39 non-vacuous steps of 8 727 admitted; 3 are refused and the option-on walk loses exactly 3 nodes, all in another SNF class (the E-332 crossing), 0 in-class nodes and 0 LNAs. By the pre-registered threshold the answer is NO (the option loses no join) for this population, on a denominator of 39 non-vacuous steps and 3 refusals; that is a small denominator, not a law. The n = 9 `merges.py` walks are VACUOUS: 0 of 7 262 distinct admitted steps has d >= 2 (even 3033030, which reaches parallel arrows). The E-316 joins are vacuous as a loss test: all 20 have both ends at the same LNA (they show the J != 0 child lies in the root's component), so refusing their audited step disconnects nothing; at path level the option refuses exactly that step in each (20 of 20) and no other. Refuted by: an in-class node or LNA lost at another E-332 root or deeper, or a refused step with d <= 1 (would contradict E-338).

## Pre-registration (written before running; not edited after)

Question: does `refuseNonTilting=True` lose joins on walks that hold steps where it can fire?

- Instrument: `workshop/rounds/133/toolsmith_dcount.py`. At every gate-admitted step of the walk it records
  d = max_{i != v} dim e_i A e_v, `isTiltingAtVertex`, legality and whether the child keeps the root's Coxeter key.
- Non-vacuous step := gate-admitted, legal, key-kept, d >= 2 (by E-338 the only steps where the option can refuse
  and change the walk). The denominator reported is the number of such steps.
- Populations: (a) library walks (`mutationSearchDepthFirst`, defaults, `fingerprint.Visited`) at E-332's two
  crossing roots (303500030, dual of 303305000), depth 6; (b) `merges.py`'s n = 9 searches (every leftover
  member and its relation dual, depth 5 as E-036); (c) recorded mixed F/R paths: E-316's joins as replayed in
  E-321 (14 + 5 regenerated + E-255) and the E-332 10-step walk (path level).
- Loss unit, library walks: an LNA reached by the default walk and not by the option-on walk from the same root at
  the same depth. Recorded paths: a join whose two ends are DISTINCT LNAs and that contains a refused step.
- In-class := same SNF(C + C^T) and signature as the root; else cross-class.
- Threshold: "the option loses joins" = YES if >= 1 in-class loss in any population. NO if 0 in-class losses and
  the population has >= 1 non-vacuous step; with 0 non-vacuous steps the population is reported VACUOUS (no answer).
  A cross-class loss is the intended effect and is reported separately.
- Control (pass is void if it fails): E-318's planted step (vertex 3 of `tests/test_search_soundness._e318Parent`)
  must have `isTiltingAtVertex` False, d = 2, and the option-on library walk of depth 1 from P must not reach its child.
- Also counted: non-vacuous steps that are tilting (d >= 2 with J = 0; the converse of E-338 fails there).

## Evidence

Control (pass valid): E-318 parent, vertex 3: gate True, tilting False, d = 2, legal, key kept; depth-1 walk from P reaches the child with the option off (8 nodes), not with it on (7 nodes).

| population | admitted steps | d >= 2 | non-vacuous | refused (kept) | refused (key-moved, vacuous) | non-vacuous tilting | lost in-class | lost cross |
|---|---|---|---|---|---|---|---|---|
| lib 11 d6 303500030 | 3 825 | 40 | 38 | 2 (d = 2, v = 8) | 2 | 36 | 0 nodes, 0 LNAs | 2 nodes |
| lib 11 d6 dual(303305000) | 4 902 | 3 | 1 | 1 (d = 2, v = 7) | 2 | 0 | 0 nodes, 0 LNAs | 1 node |
| merges n = 9 d5 (9 distinct walks) | 7 262 | 0 | 0 | 0 | 0 | 0 | 0 LNAs (VACUOUS) | 0 |
| recorded paths (22 walks, 280 steps) | 280 | 23 | 23 | 22 (20 at d = 2, 2 at d = 3) | 0 | 1 | 0 (ends = same LNA in 20 of 20 joins) | E-318, E-332 steps |

- Node-level (post hoc, not pre-registered; `nodes` mode, canonicalKey): 1 844 -> 1 842 and 2 437 -> 2 436 nodes, lost paths `1 3 4 1 3 8`, `1 1 4 3 1 8` and `11 6 11 7 7 7`, i.e. E-332's P_a / P_b crossing steps; nodes in another SNF class reached by default: 2 and 1, all lost with the option; nothing gained. The LNA-level unit I pre-registered is weak here: both depth-6 walks reach only their root as an LNA, with or without the option.
- No refused step has d <= 1 anywhere (0 of 16 269 distinct counted steps: 8 727 + 7 262 + 280), as E-338 predicts.
- d >= 2 alone does not make a step non-tilting: 37 of 62 non-vacuous steps are tilting (36 at 303500030, 1 in E-234 4551). This fits E-338, which also needs a non-monomial piece of I, so nothing new.
- The n = 9 rows come in relation-dual pairs (3505000/5050030, 3345000/5040330, 3455000/5504030 swap root/rdual with equal counts); 4550400, 4444400 and 3033030 are self-paired. 14 524 rows, 7 262 distinct.
- Cost of the option: walks with it on take as long as the default (198 vs 198 s, 206 vs 248 s, n = 9 within +-15 %).

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/133/toolsmith_dcount.py control            # 3 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/133/toolsmith_dcount.py lib 11 6 303500030  # 400 s (off counted + on)
timeout 10m .venv/bin/python -W ignore workshop/rounds/133/toolsmith_dcount.py lib 11 6 303305000 dual   # 455 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/133/toolsmith_dcount.py nodes 11 6 303500030      # ~7 min (post hoc)
timeout 10m .venv/bin/python -W ignore workshop/rounds/133/toolsmith_dcount.py lib 9 5 <member> root|rdual   # 20-45 s each; members: `members 9`
timeout 10m .venv/bin/python -W ignore workshop/rounds/133/toolsmith_dcount.py paths              # 5 s
```
All outputs: `workshop/rounds/133/toolsmith_dcount.out` (9 KB). Three jobs at a time on 4 cores, so timings are upper bounds.

## Prior record

E-332 already has "refuseNonTilting removes the crossing (2 of 2 roots)" and 3 824 -> 3 822 visits; this repeats it (my node-level loss count is that fact again). E-316 / E-321 already say each of the 19 joins has exactly one J != 0 step and the child lies in the parent's J = 0 component; my 20 of 20 single refusals match. E-036 already covers the n = 9 walks (depth 5, 3033030 the only one reaching parallel arrows) but has no d count. New: the counts of non-vacuous steps (39 / 0 / 23), 0 in-class losses at E-332's roots, n = 9 merges vacuous, and the observation that the E-316 joins cannot serve as a loss test (same-LNA ends). E-337 is still not tested on a non-thin population: this pass is the first one with d >= 2 steps in a library walk, and it is two roots.

## Code changed

None in the library. New script `workshop/rounds/133/toolsmith_dcount.py` (wraps `mutation.mutationIsPossibleAtVertex` with a counter; sympy SNF, since python-flint is not installed in `.venv`, so `maverick_libwalk.py` does not run here). Ran `tests/test_search_soundness.py -m "not slow"`: 7 passed (unchanged file, checked for the E-318 fixture).

## Next

- Chair / skeptic: with 0 in-class losses on 39 non-vacuous steps, the threshold gives NO, but the denominator is two roots. The cheap next denominator is the other 5 E-332 13-class roots at depth 6 (about 7 min each, < 10 min per command), before any n = 11 depth-4 F-041 run.
- Toolsmith: a planted in-class loss control (a refused in-class step forced by monkeypatching the test), so the loss side has power; and a mixed F/R walk with the option, since E-316/E-318 walks are not library walks.
- Theorist: the d >= 2 tilting steps (36 at 303500030) are a ready population for "d >= 2 + non-monomial piece => J != 0".
