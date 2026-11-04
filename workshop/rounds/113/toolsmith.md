# Hom-verified tilting chains join parent and child for 12 of E-201's 13 key-kept J != 0 steps at n = 7 class 1; step 8 not reached

author: toolsmith · round: 113 · kind: result
thread: T10 · bears on: H-015, R-candidate "a key-guarded walk stays in one class", E-201, E-212, E-295
scope: n = 7, key class 1, the 13 D = 0 key-kept J != 0 steps from `069/skeptic_c1_13.py` (checkpoint nexp 23000, steps indexed 0..12 in that order). Bidirectional search of vertex tilting edges (F/R at one vertex), depth <= 5 per side (steps 3, 9: 5 and 6), cap 40000 (steps 3, 9: 150000) algebras; search uses the repo's gate/J/key filters only to propose edges, verification uses none. Generation of K^b(proj) by T is only proxied by the K0 determinant (Okuyama-Rickard assumed, as in E-168). Meets are structure-graph isomorphisms (quiver, arrow multiplicities, relation shapes), relation coefficients are not compared.
power: verifier says NO on 3 planted classes: (P1/P3 by the skeptic's design) corrupted single edges and wrong end; here, per step, 78-104 single-edge corruptions (steps 2, 4, 5): all rejected except 2 each that are themselves valid tilting edges into an isomorphic algebra (true YES, not false positives); 12 of 12 children of other steps substituted for the end: rejected (36 of 36 over steps 2, 4, 5); P2 duplicated summand: K0 det 0 in 708 of 708 (real T: det +-1 in 118 of 118). It has NO power for "child leaves the class": no key-equal inequivalent pair at n <= 7. An absent chain is UNDECIDED.

## Claim
For 12 of the 13 steps (all but 8) there is an explicit chain parent = X_0 - ... - X_a ~ X_b - ... - X_0' = child of 6 to 11 vertex tilting edges, every edge having Hom_K(T,T[-1]) = Hom_K(T,T[1]) = 0 (chain maps mod homotopy, `050/skeptic_tilt.py`), K0 determinant of the 7 summands +-1, and dim Hom(T_i,T_j) equal to the labelled Cartan matrix of the next algebra, the two halves meeting in an isomorphic presentation. So parent and child are derived equivalent for those 12 (modulo generation, assumed). The direct step parent -> child itself FAILS the Hom test in all 13 (dim Hom(T,T[-1]) = 1; this is E-201, expected, not a NO for equivalence, E-212). Step 8: no verdict (see below). The claim "13 of 13" is NOT made; the R-candidate stays open on one step (8). A NO would have been: a chain edge with Hom(T,T[+-1]) != 0 or K0 det != +-1 or Cartan mismatch, or a non-isomorphic meet; none occurred in the 12 chains (all edges passed). Not a NO for any step: an exhausted search is only a bound.

## Evidence
| step | v | parent/child parallel arrows | edges (A-side | B-side) | verdict |
|---|---|---|---|---|
| 0 | 4 | yes/yes | F1 F3 F7 R5 R4 | R4 F7 F1 F3 | YES (100 s) |
| 1 | 4 | yes/yes | F5 F7 F1 R3 R4 | R4 F1 F5 F7 | YES (196 s) |
| 2 | 6 | no/no | R7 R6 R3 | F4 R6 F1 | YES (44 s) |
| 3 | 7 | no/no | R2 F3 F2 F3 F7 | F3 F5 R3 F6 R5 F3 | YES (11 edges; needed depth 5 + 6; 151 s on resumed ball 17 k / 10 k) |
| 4 | 4 | no/no | F3 R5 R4 R3 | F3 F7 R4 F1 | YES (73 s) |
| 5 | 4 | no/no | F1 R5 R4 R3 | F1 F7 R4 F1 | YES (69 s) |
| 6 | 4 | yes/yes | F1 F3 F7 R5 R4 | R4 F7 F1 F3 | YES (120 s) |
| 7 | 5 | yes/yes | F1 F3 F7 R4 R5 | R5 F7 F1 F3 | YES (107 s) |
| 8 | 4 | yes/yes | | | not reached: no meet; ball 5.5 k (A) / 9.8 k (B) algebras after about 4 x 8 min, parallel-arrow algebras are slow; UNDECIDED |
| 9 | 7 | no/no | R2 F1 F2 F1 F7 | F1 F5 R1 F6 R5 F1 | YES (11 edges; depth 5 + 6; 180 s on resumed ball 21 k / 7.5 k) |
| 10 | 1 | no/no | F6 F7 R3 F6 F3 | F6 F6 R1 R6 F7 | YES (290 s) |
| 11 | 1 | no/no | F6 F4 R7 F6 F7 | F4 F6 F6 R1 R6 | YES (280 s) |
| 12 | 7 | no/no | R1 F4 F1 F4 F7 | F2 F4 F2 F4 R2 | YES (276 s) |
Edge letters: F = right-mutation-as-coded on the algebra, R = on the dual; vertex number. Steps 0 and 6 got the same move strings (1 and 7 are the same up to relabelling). All meets are in algebras without parallel arrows, so the isomorphism is an honest one up to relation coefficients. Steps 3 and 9 (same vertex 7) are mirror images of each other (R2 F3 F2 F3 F7 vs R2 F1 F2 F1 F7) and need one more level on the B side.
Controls in Power. Corruption "accepted 2 of 104": the replaced last edge is itself a valid tilting edge to an isomorphic algebra, so those are true YES.

## Reproduction
```
timeout 10m .venv/bin/python -u workshop/rounds/069/skeptic_c1_13.py 1 540 23000    # twice; makes /tmp/skeptic_c1_13_ck1.pkl, ~9 min each
timeout 10m .venv/bin/python -u workshop/rounds/113/toolsmith_c1_chain.py list
timeout 10m .venv/bin/python -u workshop/rounds/113/toolsmith_c1_chain.py chain I 5 40000   # I = 0..12 (steps 3, 8, 9: LIM=5,6 ... chain I 6 150000); rerun if it prints PAUSED (checkpoint /tmp/c1chain_ckI.pkl)
timeout 10m .venv/bin/python -u workshop/rounds/113/toolsmith_c1_chain.py ctl 4 4 12000     # corruption and end-swap controls
timeout 10m .venv/bin/python workshop/rounds/113/toolsmith_c1_k0ctl.py                      # K0 control, about 2 min
```

## Prior record
E-201 (the 13 steps, Hom(T,T[-1]) != 0), E-295 (Cartan route cannot decide, 89 steps congruent), E-212 (class 2: 8 of 9 Hom-verified chains, record 6 open), E-161/E-163 (Hom test), E-168 (generation per step). New: chains for 12 of the 13 n = 7 class-1 steps with the K0 check added to the verifier (E-212's verifier had none) and a K0 power control.

## Code changed
`workshop/rounds/113/toolsmith_c1_chain.py` (finished the earlier partial: added `k0det`, folded into `horn`), `toolsmith_c1_k0ctl.py`. No library changes, no tests needed.

## Next
Step 8: a deeper search (resumable; proposal for OVERNIGHT.md: `LIM=6,7 .venv/bin/python -u workshop/rounds/113/toolsmith_c1_chain.py chain 8 7 400000` repeated, each call 8 min) or a cone-closure approach as `107/toolsmith_cones.py`. Skeptic: whether 12 of 13 plus 1 undecided is enough to move the R-candidate.
