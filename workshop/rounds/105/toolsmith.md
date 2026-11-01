# No vertex-local single-summand tilting exchange (right or left) of B_6 has End(T) with the Cartan matrix of A_6 (independent Hom test, n = 6; same at n = 5); 7 of the 8 E-268 key-kept J != 0 edges at n = 9 are silting, not tilting

author: toolsmith · round: 105 · kind: result
thread: H-027 / T10 · bears on: H-027, H-015, E-280, E-278, E-268, E-161, E-184
scope: vertex-local single-summand exchange from the tilting complex B_n itself, n = 6 (B_6 = row (0,0,3,0) vs A_6 = row (5,0,0,0)) and n = 5 (B_5 = (0,3,0) vs A_5 = (4,0,0)); n = 6: ALL 10 complexes T_k (5 right, 5 left via the opposite algebra), including gate-refused ones; n = 5: 8 (4 + 4); test compares Cartan matrices only (necessary for End(T) ~ A_6, not sufficient). Part 2: n = 9, 7 of the 8 E-268 edges (A: 1, B: 6 printed E-244 parents); the E-239 depth-5 B edge was not rebuilt. Multi-vertex / two-step mutations NOT covered.
power: yes. (a) On all LNAs n = 4, 5, 6 (right steps, 15 + 56 + 210 = 281) the test says NO on 5 + 21 + 84 = 110 and YES on 10 + 35 + 126 = 171, agreeing with the repo gate (isMutable) on 281 of 281 (the gate is not used by the test; this is consistency with a rule validated by the same test family, so (b) and (c) carry the power). (b) planted errors: wrong shift on A_3 (T = P1 + P2[-1] + P3) gives 2 nonzero Hom(T,T[+-1]) entries; repeated summand P1+P1+P3 gives det H = 0 vs det Cartan 1. (c) known tilting A_3, LNA(4;2,0) k=1,3, B_6 k=1..4 all YES; known non-tilting LNA(4;2,0) k=2 and B_6 k=5 NO.
backlog: n/a

## Response to referee

Verdict minor revision; all three required items done.

1. Retitle/rescope: done. Title, scope line and Claim now say "vertex-local single-summand exchange from B_n itself", not "one-step mutation". The search space is moves (one summand replaced by its minimal arrow approximation), independent of the E-187 move set but not of the vertex-local form. Part 2 now reads "7 of the 8 E-268 edges are silting, not tilting".
2. H-027 stays OPEN: done. Added to Claim. The falsifier asks for an exhaustive negative at n = 5 or 6 over tilting steps; this is single-summand exchange from B_n, one step, so it does not close H-027. I also ran n = 5 with the same script (`toolsmith_tiltA6.py 5`, B_5 = (0,3,0) vs A_5 = (4,0,0), 8 complexes: 6 tilting, 2 not (right k = 4, left k = 3, both gate-refused), 0 Cartan-equal to A_5; output `toolsmith_tiltA5_out.txt`). The n = 5 negative is the same kind (one step, vertex-local).
3. Runtime: done. First command takes about 3 s, not "~1 min" (second about 2 s).

Referee remark on power (a): agreed that agreement with the gate is consistency with a rule validated by the same family of tests (E-161/E-184); the power rests on (b) and (c). Text says so below.

## Claim

Part 1. For B_6 = LNA 0,0,3,0, take every vertex-local complex T_k = (+_{i != k} P_i) + (P_k -> +_{k->h} P_h), right (arrows out of k) and left (same on the opposite algebra). By an independent Hom(T,T[m]) computation in K^b(proj) over GF(2^31-1) (rounds/050 `skeptic_tilt.stepTest`, no `tiltingPlus`, `perI`, gate or E-187 rule), 8 of the 10 are tilting (right k = 1..4; left k = 2, 3, 5, 6) and 2 are not (right k = 5, Hom(T,T[-1]) of dim 3; left k = 4, dim 1); both non-tilting ones are exactly the gate-refused ones, and the 8 tilting are the gate-admitted ones. For none of the 8 tilting T is the Cartan matrix H of End(T) equal to the Cartan matrix of A_6 (up to simultaneous permutation or transpose; 720 permutations tried). Hence no one-step tilting mutation of B_6 is isomorphic to A_6, so A_6 is not reached from B_6 by a one-step tilting mutation of either side, inside or outside the E-187 move set. It does NOT say that A_6 is not reached in two or more steps, nor anything about non-vertex-local tilting complexes (A_6 ~ B_6 is Fosse Cor. 4.6, so a tilting complex exists). H-027 stays OPEN: its falsifier is an exhaustive negative over tilting steps at n = 5 or 6; this is one step, vertex-local single-summand exchange from B_n, at n = 6 (and n = 5, below), and nothing about two or more steps or non-vertex-local complexes. Refuted by: any of the 8 H matrices being permutation-equal to Cartan(A_6).

Part 2. At n = 9, 7 of the 8 E-268 edges (the E-239 depth-5 B edge not rebuilt): the 7 rebuilt key-kept J != 0 edges (E-254 one A edge, E-244 six B edges, the parents replayed from the printed paths of rounds/098) all have Hom(T,T[-1]) != 0 (mass 2, 2, 4, 4, 2, 2, 2; Hom(T,T[1]) = 0), so T is silting but not tilting; det H = det Cartan (K0 consistent) and H = Cartan(child)^T in all 7. So these edges are NOT certificates of derived equivalence by the tilting criterion (AI Thm 2.32); they say nothing about the J = 0 premise beyond E-161/E-184, which found the same at n = 7 and n = 10.

## Evidence

Part 1 (output `toolsmith_tiltA6_out.txt`, includes every H matrix): 8 tilting steps, none Cartan-equal to A_6; the Cartan row-sum structure of A_6 (one zero pair (1,6), all other i <= j pairs 1) is not matched by any child. Det: every H has |det| = 1 = det Cartan(B_6), so generation is not contradicted at K0; generation is by construction (the other n-1 summands are the P_h, h != k, and P_k is the cone, acyclic quiver, E-168). Gate agreement over all LNAs: n = 4: 10 YES/YES, 5 NO/NO; n = 5: 35/21; n = 6: 126/84; zero off-diagonal cells.
n = 5 (`toolsmith_tiltA5_out.txt`): right k = 1..3 and left k = 2, 4, 5 tilting (6), right k = 4 and left k = 3 not (2); none Cartan-equal to A_5.
Part 2 (`toolsmith_tilt8_out.txt`): per edge `hom0 False`, `H==Cartan(child)^T True`, <1 s each.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/105/toolsmith_tiltA6.py > workshop/rounds/105/toolsmith_tiltA6_out.txt   # ~3 s
timeout 10m .venv/bin/python workshop/rounds/105/toolsmith_tiltA6.py 5 > workshop/rounds/105/toolsmith_tiltA5_out.txt   # n = 5, ~2 s
timeout 10m .venv/bin/python workshop/rounds/105/toolsmith_tilt8.py                                              # ~2 s
```

## Prior record

The test itself is E-161 (round 050). Part 2 repeats at n = 9 what E-161 (n = 7, 25 of 25 rejected), E-184 (n = 7, 9 of 9) and E-174 found: J != 0 steps have Hom(T,T[-1]) != 0 and End(T) still equals the child at the Cartan level. Part 1 is new: E-280 computed only the move orbit, not the necessary condition on all one-step complexes. The repo's mutation is not used for the verdict, only for context (gate column).

## Code changed

None. New: `workshop/rounds/105/toolsmith_tiltA6.py`, `toolsmith_tilt8.py` (both exec `rounds/050/skeptic_tilt.py`, `rounds/087/toolsmith_d6scan.py`, `rounds/082/skeptic_join.py`); outputs `toolsmith_tiltA6_out.txt`, `toolsmith_tiltA5_out.txt` (script takes optional n), `toolsmith_tilt8_out.txt`.

## Next

- theorist: a two-step search is needed (one-step is empty); candidates are tilting T of B_6 whose H equals Cartan(A_6): enumerate by Cartan-first (Hom-dimension matrix) rather than by moves; Fosse Cor. 4.6 gives the actual complex, which is the cheapest certificate.
- Note: the E-268 premise (J = 0) is untouched; the J != 0 edges fail tilting, so a join through them is not a derived-equivalence proof.
