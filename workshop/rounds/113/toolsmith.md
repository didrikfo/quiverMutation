# Hom-verified tilting chains, with meets isomorphic up to relation coefficients and generation assumed (E-168), join parent and child for 12 of E-201's 13 class-1 steps at n = 7; a YES-side consistency check with no NO case available; step 8 not reached

author: toolsmith · round: 113 · kind: result
thread: T10 · bears on: H-015, R-candidate "a key-guarded walk stays in one class", E-201, E-212, E-295
scope: n = 7, key class 1, the 13 D = 0 key-kept J != 0 steps from `069/skeptic_c1_13.py` (checkpoint nexp 23000, steps indexed 0..12 in that order). Bidirectional search of vertex tilting edges (F/R at one vertex), depth <= 5 per side (steps 3, 9: 5 and 6), cap 40000 (steps 3, 9: 150000) algebras; search uses the repo's gate/J/key filters only to propose edges, verification uses none. Generation of K^b(proj) by T is only proxied by the K0 determinant (Okuyama-Rickard assumed, as in E-168). Meets are structure-graph isomorphisms (quiver, arrow multiplicities, relation shapes), relation coefficients are NOT compared (a sign or scalar in a two-term relation could change the algebra). Every mutated vertex of the 12 chains is loop-free (checked in round 113 response, see below).
power: verifier says NO on 3 planted classes, split by the test that rejected (steps 2, 4, 5; 78 + 104 + 104 = 286 single-edge corruptions): accepted 6 (2 per step; each is itself a valid tilting edge into an isomorphic algebra, a true YES, not a false positive); rejected 280, of which meet iso rejected 280 (every one: a failed step makes the meet empty), Hom(T,T[+-1]) != 0 additionally in 11 (3+4+4), "step impossible / reduction changes vertex set" in 151, K0 det 0 times, Cartan 0 times, loop 0 times. Only the meet iso rejected 35+40+45 = 120 corruptions that no other test flagged. So (as in E-228) Hom and Cartan alone have almost no NO power on corruptions; K0 and Cartan have none here. 12 of 12 children of other steps substituted for the end (36 of 36 over steps 2, 4, 5): all rejected by the meet iso alone (some also by the empty-step test). P2 duplicated summand is caught by K0 only: det 0 in 708 of 708 (real T: det +-1 in 118 of 118; separate object set, `toolsmith_c1_k0ctl.py`). Not tested: controls for the other 9 chains (not run). There is NO power for "child leaves the class": no key-equal inequivalent pair at n <= 7. An absent chain is UNDECIDED.

## Response to referee
Verdict received: minor revision. Items 1-4 done; item 5 (step 8 deeper search) is next round, not done.
Re-runs (this round, 4 cores, 9 jobs at once so the times in the table are not comparable): chain re-run with the new loop test for 0, 1, 2, 4, 5, 6, 7, 10, 11, 12: all YES, same move strings as the table. Steps 3 and 9: the resumed-checkpoint search returned DIFFERENT 11-edge chains (3: F5 F3 R2 F5 F2 F5 | F5 F6 F5 F3 R7; 9: F5 F1 R2 F5 F2 F5 | F5 F6 F5 F1 R7), verdict YES but with parallel arrows at the meet, where the iso is only a necessary test. The table's original chains for 3 and 9 were therefore replayed with the new verifier (`chain`-file mode `replay`): both YES, meet without parallel arrows. The table is unchanged.
1. Claim, title and Scope now say: meet iso ignores relation coefficients, generation by T assumed (E-168), derived equivalence "up to coefficients".
2. Loop-freedom: the verifier now tests number_of_edges(v,v) == 0 at every mutated vertex of every edge (flag `!` in the log). Result: loop-free at all edges of all 12 chains (0 flags; 3 and 9 via replay of the original chains plus the re-found ones). So E-168 generation applies at each edge; the K0 determinant is then a cross-check, not the only argument.
3. Power split by test: see the power line. The meet iso carries all rejections, Hom 11 of 280, K0/Cartan/loop 0. A YES here is carried by Hom(T,T[+-1]) = 0 on every edge (which is a real tilting statement given generation, and every edge passed it) plus the meet iso.
4. This is a YES-side consistency check: no NO case exists at n <= 7 (no key-equal inequivalent pair), so it could not have refuted the R-candidate. It does not move the R-candidate beyond E-228's level (12 of 13 here, no NO power). Said in the Claim.
5. Step 8: next round (the overnight proposal under Next stands).

## Claim
For 12 of the 13 steps (all but 8) there is an explicit chain parent = X_0 - ... - X_a ~ X_b - ... - X_0' = child of 6 to 11 vertex tilting edges, every edge having Hom_K(T,T[-1]) = Hom_K(T,T[1]) = 0 (chain maps mod homotopy, `050/skeptic_tilt.py`), loop-free mutated vertex, K0 determinant of the 7 summands +-1, and dim Hom(T_i,T_j) equal to the labelled Cartan matrix of the next algebra, the two halves meeting in an isomorphic presentation. So parent and child are derived equivalent for those 12, up to relation coefficients (meet iso ignores them) and assuming generation by T (E-168: every mutated vertex loop-free, verified; Okuyama-Rickard). This is a YES-side consistency check: no NO case is available at n <= 7 (no key-equal inequivalent pair), so it could not have refuted the R-candidate and does not move it beyond E-228's level (Hom + Cartan had no NO power there, 0 of 52; here the meet iso did all the rejecting, Hom 11 of 280). The direct step parent -> child itself FAILS the Hom test in all 13 (dim Hom(T,T[-1]) = 1; this is E-201, expected, not a NO for equivalence, E-212). Step 8: no verdict (see below). The claim "13 of 13" is NOT made; the R-candidate stays open on one step (8). A NO would have been: a chain edge with Hom(T,T[+-1]) != 0 or K0 det != +-1 or Cartan mismatch, or a non-isomorphic meet; none occurred in the 12 chains (all edges passed). Not a NO for any step: an exhausted search is only a bound.

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
timeout 10m .venv/bin/python -u workshop/rounds/113/toolsmith_c1_chain.py ctl 4 5 40000     # corruption and end-swap controls, split by test
timeout 10m .venv/bin/python -u workshop/rounds/113/toolsmith_c1_chain.py replay 3 0 0 "R2 F3 F2 F3 F7" "F3 F5 R3 F6 R5 F3"   # explicit chains of steps 3 and 9 (table)
timeout 10m .venv/bin/python workshop/rounds/113/toolsmith_c1_k0ctl.py                      # K0 control, about 2 min
```

## Prior record
E-201 (the 13 steps, Hom(T,T[-1]) != 0), E-295 (Cartan route cannot decide, 89 steps congruent), E-212 (class 2: 8 of 9 Hom-verified chains, record 6 open), E-161/E-163 (Hom test), E-168 (generation per step). New: chains for 12 of the 13 n = 7 class-1 steps with the K0 check added to the verifier (E-212's verifier had none) and a K0 power control.

## Code changed
`workshop/rounds/113/toolsmith_c1_chain.py` (finished the earlier partial: added `k0det`, folded into `horn`), `toolsmith_c1_k0ctl.py`. No library changes, no tests needed.

## Next
Step 8: a deeper search (resumable; proposal for OVERNIGHT.md: `LIM=6,7 .venv/bin/python -u workshop/rounds/113/toolsmith_c1_chain.py chain 8 7 400000` repeated, each call 8 min) or a cone-closure approach as `107/toolsmith_cones.py`. Skeptic: whether 12 of 13 plus 1 undecided is enough to move the R-candidate.
