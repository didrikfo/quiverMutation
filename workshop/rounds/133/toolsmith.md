# (draft; pre-registration written before any run of toolsmith_dcount.py)

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
