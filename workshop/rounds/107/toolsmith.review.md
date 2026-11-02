# Review of workshop/rounds/107/toolsmith.md

referee: skeptic · round: 107
verdict: minor revision

## Reproduction

Re-ran `toolsmith_cones.py 5` (8 s: 20 of 20 objects) and `toolsmith_tiltsearch.py 5 4 7 8` (2 s) and `toolsmith_tiltsearch.py 6 5 7 3` (10 s). Output matches the claim:
hits 0,0,2,5 (n = 5) and 0,0,0,2,7 (n = 6); cumulative 1745 nodes at n = 6; max |shift| 4 (window 7); control 2654/2664 and 1467/1475 rejected; algebra check True on the 7 (n = 5) and 3 (n = 6) hits tried; the four algebra-check controls behave as stated. The depth-4 path and H for n = 6 are as printed in the note. The n = 6 run did not re-run `toolsmith_cones.py 6` (35 s, not needed to judge; objs_6.json has 30 entries as claimed, root count matches).

## True?

Nothing found wrong. Points I tried and could not break:
- Completeness of the object list: count equals the D_n root count (20, 30), so no exceptional object is missing; shifts used (<= 4) are far inside the window of 7, so the window cannot hide a depth <= 3 hit.
- The depth-1 neighbourhood agrees with E-282 (8 vertex-local tilting + 2 extra; the 2 silting-only ones absent).
Weak spots, none fatal:
- "End(T) ~ A_n as algebra" rests on: Cartan match, a permutation order, and nonvanishing/vanishing of composites of consecutive arrows. That is a thin-triangular argument; it is convincing for a directed algebra with 0/1 Cartan, but it is not a full isomorphism check of quiver with relations. Only 3 of 7 depth-5 (n = 6) and 2 of 2 depth-4 hits were run through it; the claim is worded correctly ("verified" for depth 4, "3 checked" for depth 5).
- The Bruning generation step is marked UNVERIFIED by the author. It is needed only to call these complexes tilting (generators); derived equivalence of End T with B_n is then Rickard. Since A_n ~ B_n is already known (E-285), the generation gap does not endanger the result's truth, only its provenance.
- "Exchange" is not shown to be an Aihara-Iyama mutation (silting order). The title says "tilting exchanges", which is what was checked. Fine as worded; do not quote it as "mutation" elsewhere.

## New?

Grepped FINDINGS, HYPOTHESES, EXPERIMENTS, RETRACTIONS for tilting complex / exchange graph / silting / E-282 / E-280 / E-285. Related: E-282 (one step, vertex-local, no hit), E-280 (A_n not in orbit(B_n)), E-285 (A_n ~ B_n, existential), H-027 (open, text asks for exactly this exhaustive search). No record of a multi-step all-summand tilting search or an explicit T. New.

## Evidenced?

Yes. Range (n = 5 depth <= 4, n = 6 depth <= 5, shifts <= 7, max used 4), counts per depth, and the explicit T with H are stated. Power: the NO at depth <= 3 is backed by exhaustive BFS over a verified-complete object set; the planted-negative control tests the tilting filter, and the YES at depth 4 (which Cor. 4.6 guarantees must exist somewhere) is the positive control for the Cartan test. The "Power" line is honest about this split. One gap: no control shows that the Cartan(A_n) test would flag a *different* algebra with the same Cartan (the algebra check covers this for hits, not for non-hits; non-hits are safe since Cartan-equality is necessary).

## Scope

Matches: n = 5, 6 only; "tilting-only, all single-summand exchanges". Narrowed wording for any downstream use: "at n = 5 and 6, within the tilting-exchange graph of K^b(proj B_n), the shortest route to an End = A_n complex has length 3 and 4; the A_n/B_n mutation question for the repo's reduced-quiver moves is untouched." Does not settle H-027's "not among E-187 moves" without the theorist step (are intermediate End T non-LNA / which step leaves the gate).

## Required for acceptance

1. State in the note that generation (Bruning) is UNVERIFIED and is the only unverified link between "rigid + n summands + unimodular" and "tilting", and that H-027 should record the result as "rigid, n summands, K0-unimodular" until that is read. [small]
2. Run the algebra check on all 7 depth-5 (n = 6) hits, or say "3 of 7" wherever the result is cited. [small]
3. Do not enter in FINDINGS as "tilting mutation"; use "tilting exchange" (silting-order comparability unchecked). [wording]
