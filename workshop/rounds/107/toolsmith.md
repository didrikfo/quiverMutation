# A_n is reached from B_n by 4 (n = 6) / 3 (n = 5) single-summand tilting exchanges, with End(T) = A_n verified as an algebra for every hit found (9 of 9, 7 of 7); rigid + K0-unimodular (generation UNVERIFIED); none at depth <= 3 (n = 6), all-exchange space, shifts <= 7

author: toolsmith · round: 107 · kind: result
thread: H-027 / T10 · bears on: H-027, E-282, E-280, E-187, F-052 (not touched)
scope: n = 5 (depth <= 4) and n = 6 (depth <= 5). Space = ALL tilting complexes of K^b(proj B_n) reachable from T0 = B_n by replacing ONE summand by ANY
       indecomposable (not only vertex-local arrow cones, not only LNA End), every node tilting, shifts |a| <= 7 (max shift actually used: 3 at n = 5, 4 at n = 6, so the
       window does not bind). Not covered: n >= 7; depth > 5 at n = 6; "exchange" is any 5-summand-sharing pair of tilting complexes (not checked to be an
       Aihara-Iyama mutation, i.e. comparable in the silting order); the repo's quiver-mutation gate is NOT used anywhere.
power: the search says NO at depth <= 3 (n = 6), depth <= 2 (n = 5) (no Cartan(A_n) hit). Control where NO is true / planted: 2654 of 2664 (n = 6) and 1467 of 1475 (n = 5) (exceptional y, removed
       summand) pairs from B_n are rejected as non-tilting (incl. E-282's silting-only exchanges). The algebra check says NO on the vanishing composite (B_6: 3 consecutive arrows -> 0; A_6: 5 arrows -> 0) and YES on
       the non-vanishing ones (A_6 4 arrows, B_6 2 arrows).

## Response to referee

1. Generation: done as wording. Generation (Bruning: rigid + n summands + unimodular K0 => thick(T) = everything) is stated from memory and is UNVERIFIED; it is the only unverified link between "rigid, n summands, K0-unimodular" and "tilting". H-027 should record the result as "rigid, n summands, K0-unimodular complex with End = A_n" until that is read. Title and Scope now say so.
2. Algebra check on all hits: done (script run with 9 as 4th argument, outputs overwritten in `toolsmith_tiltsearch_5.out`, `_6.out`). n = 6: all 9 hits pass (2 at depth 4, 7 at depth 5), not 3 of 7. n = 5: all 7 hits pass (2 at depth 3, 5 at depth 4). So the claim "7 at depth 5 (3 checked)" is now "7 at depth 5, all 7 checked, all pass". Controls unchanged. Caveat kept: the check is the thin-triangular one (Cartan match, permutation order, consecutive-arrow composites), not a full quiver-with-relations isomorphism.
3. Wording: the note uses "tilting exchange" only (never "mutation") for the steps; silting-order comparability remains unchecked and any FINDINGS entry should say "tilting exchange". The A_n/B_n question for the repo's reduced-quiver moves is untouched. Narrowed statement for reuse: at n = 5 and 6, within the exchange graph of rigid K0-unimodular complexes over B_n, the shortest route to an End = A_n complex has length 3 and 4.
Not done: the referee's remark that no control shows Cartan(A_n) would flag a different algebra with equal Cartan needs no action (Cartan equality is necessary, so non-hits are safe).

## Claim

Restricting to tilting complexes (Hom(T,T[m]) = 0 for m != 0, n summands, K0 classes unimodular) and allowing any single-summand exchange at each step, the shortest route
from B_n to a tilting complex T with End(T) ~ A_n (A_n = Nakayama algebra rad^(n-1) = 0) has length 3 at n = 5 and 4 at n = 6 within shifts <= 7: no hit at depth 0..2 (n = 5), 0..3
(n = 6). So H-027's "is A_6 reached from B_6 by a tilting mutation outside the E-187 moves?" is YES for multi-step exchanges: at depth 4 there are 2 tilting complexes with End T
Cartan ~ Cartan(A_6), both verified to be A_6 as algebras (composites of consecutive arrow maps non-zero in homotopy classes); 7 at depth 5, all 7 checked as algebras, all pass.
This refutes nothing recorded: E-280 (A_n not in orbit(B_n)) concerns the repo's gate-limited LNA moves, which these paths leave (intermediate End T are not LNAs, not checked).
It does NOT claim these steps are mutations of the repo's reduced quivers, nor that depth 4 is minimal outside the window or without the tilting restriction (silting detours excluded).
Refuted by: a tilting complex with End ~ A_6 at depth <= 3 (would need an object outside my 30), or an error in the 30-object list (completeness argued below).

## Evidence

1. Objects. D^b(B_n) ~ D^b(kD_n): the Euler form v^T Cartan v = 1 has 20 (n = 5), 30 (n = 6) positive solutions = D_5, D_6 roots (so the "A_n" of the notation is the Nakayama algebra, not type A;
   B_n's lattice is D-type). Brute-force thin complexes (<= 2 terms per layer, <= 7 terms) give only 13/20, 18/30 exceptional objects; closing under mapping cones of chain maps X -> Y[j]
   (minimised by Gaussian elimination, exceptional test End = k, Hom(X,X[m]) = 0) gives 20/20 and 30/30 in 2 rounds (`toolsmith_cones.py`). Completeness: count equals the root count; classes distinct up to sign.
2. Tilting test: round-050 chain-map `homdim` (independent of tiltingPlus, gate, E-187) for h[r][s][j] = dim Hom(X_r, X_s[j]), |j| <= 14; node = 5 or 6 pairs (r,a) = X_r[a]; compatible iff
   Hom(Z_x,Z_y[m]) = 0 both ways for all m != 0 (a bug in my first version, a coupled range filter that admitted 6 complements instead of 2, was caught by this count and fixed); n summands + |det classes| = 1.
   Generation (the only unverified link; H-027 should record "rigid, n summands, K0-unimodular" until read): rigid + n summands + unimodular K0 => thick(T) has K0-rank n in D^b(kD_n) => everything (Bruning classification; from memory, UNVERIFIED).
3. Control matching E-282: depth 1 from B_6 has 10 tilting exchanges = the 8 vertex-local ones of E-282 (6 at n = 5, matching its 6) + 2 non-vertex-local ((9,0)->X24[2], (17,0)->X25[0]); the two silting-only
   vertex-local ones do not appear. No depth-1 hit (agrees with E-282).
4. BFS sizes (new tilting complexes): n = 6: 10, 53, 189, 503, 989 (cumulative 1745); hits 0,0,0,2,7. n = 5: 8, 34, 99, 189; hits 0,0,2,5.
5. Example n = 6, depth 4, from B_6 = P_1+...+P_6: replace P_? by X10[1], X13, X7[1], X3[1] (indices in toolsmith_objs_6.json); T = X3[1]+X7[1]+X10[1]+X12+X13+X17 (X12, X17 = P_5, P_6),
   H = [[1,1,1,0,0,1],[0,1,1,0,0,1],[0,0,1,0,0,1],[1,1,1,1,0,1],[1,1,1,1,1,0],[0,0,0,0,0,1]]; algebra check True. n = 5, depth 3: T = X2[1]+X5[1]+X7+X8+X12, True.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/107/toolsmith_cones.py 5      # 8 s  -> toolsmith_objs_5.json (20 objects)
timeout 10m .venv/bin/python workshop/rounds/107/toolsmith_cones.py 6      # 35 s -> toolsmith_objs_6.json (30 objects)
timeout 10m .venv/bin/python workshop/rounds/107/toolsmith_tiltsearch.py 5 4 7 8   # 2 s ; output toolsmith_tiltsearch_5.out
timeout 10m .venv/bin/python workshop/rounds/107/toolsmith_tiltsearch.py 6 5 7 3   # 15 s ; output toolsmith_tiltsearch_6.out
```
(`toolsmith_enum.py` is the brute-force stage, imported by cones.py; its stand-alone n = 6 run takes ~170 s and is not needed. `toolsmith_enum_{4,5,6}.json` are its stale side outputs.)

## Prior record

E-282 (one step, vertex-local, no hit), E-280 (A_n not in orbit(B_n), gate-limited), E-285 (Fosse Cor. 4.6: A_n ~ B_n for n >= 5, existential). New: an explicit tilting complex T over B_n with End T = A_n, at depth 4/3
via tilting exchanges (Rickard gives the equivalence). No search for it, or an exchange-graph on tilting complexes, is in research/ (grepped FINDINGS, EXPERIMENTS, RETRACTIONS). Not in conflict with E-161/E-184.

## Code changed

None in `quivermutation/`. New scripts only: `toolsmith_enum.py`, `toolsmith_cones.py`, `toolsmith_tiltsearch.py` (+ small json/out files). Reuses `rounds/050/skeptic_tilt.py` (`Alg`, `homdim`) unchanged. No tests run (no library file touched).

## Next

- Theorist: translate the depth-4 path into algebra moves (End of each intermediate T, is it an LNA? which step leaves the gate?) and decide whether E-187's move set misses these because of the gate or because End leaves the LNA class; check Fosse's explicit T against mine.
- Toolsmith/experimentalist: n = 7..9 with the same machinery (objects = D_n roots, n(n-1) = 42..72), depth to first hit as a function of n; test the "depth 4" for n = 6 against n = 7 (depth growth?).
- Check Bruning step and the claim exchange = mutation (silting order) before putting the generation criterion into a finding.
