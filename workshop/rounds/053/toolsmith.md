# At n = 7, class 1, End(T) of the tilting complex is isomorphic (label-preserving, over the algebraic closure) to the next algebra on all 13 E-161 edges; the same holds at the J != 0 failing steps, so the End(T) comparison is blind to the J = 0 premise

author: toolsmith · round: 053 · kind: result (with a negative half)
thread: T10 · bears on: E-159, E-161, H-015 (indirectly)
scope: n = 7, class 1; the 13 edges of the E-161 path (child 14: `F1 F3 F1 F5 F4 R7 R1`, LNA/dual #9: `F7 R2 R1 F2 F7 R3`); S = -1 convention of E-159; generation by T assumed (Okuyama-Rickard); negative control on the 16 failing c1 steps (8 of them with parallel arrows, undecided); class 2 and the 3 E-158 paths not run.

## Claim
For each of the 13 edges x -> y of the E-161 path, let a = x (move F) or dual(x) (move R), v the vertex, T = P_v replaced by (P_v -> sum of P_h over arrows v -> h) in degrees (-1, 0), all other summands P_i. I computed End_K(T) = Hom_{K^b(proj)}(T, T) as an algebra over Q (basis, composition, radical filtration, quiver, relations; code independent of `mutation`, `reduction`, `tiltingPlus`) and compared it with c = `reducePathAlgebra(quiverMutationAtVertex(a, v))`, the next algebra on the path (y = c or dual(c)). **13 of 13 agree**: (1) same labelled quiver (8, 7, 8, 7, 8, 8, 8 arrows on the child side; 7, 7, 9, 9, 7, 8 on the LNA side; no parallel arrows on this path), the End(T) arrow i -> j corresponds to the arrow i -> j of c (no opposite); (2) dim Hom(T_i, T_j) = dim of paths i ~> j mod I_c for all i, j (and equals the skeptic's Hom-test matrix H of E-159); (3) there are nonzero scalars lambda_a and radical-square corrections for the arrows of c such that every relation of c vanishes in End(T) (Groebner basis of the polynomial system with z * prod(lambda) = 1 is not [1]); then K Q_c / I_c -> End(T) is a surjection between equal finite dimensions, hence an isomorphism. 9 edges are settled by arrow scalars alone; 4 (F4, R7, R1 on the child side, R3 on the LNA side) need the radical-square correction (the chosen lifts do not make the relation images proportional). Also Hom(T, T[-1]) = Hom(T, T[1]) = 0 on all 13 (skeptic's code reused, E-159), so T is a presilting complex with End(T) = the next algebra.

**What this does not show, and the negative half.** The same comparison at the 16 failing J != 0 key-keeping steps (parent -> child, c1) gives iso in the 8 decided cases and "undecided, parallel arrows" in 8. So End(T) of the 2-term complex is the mutation algebra whether or not T is tilting; the quiver-level comparison therefore cannot test the J = 0 premise, and agreement on the 13 edges adds nothing beyond Hom(T,T[m]) = 0 (E-159). It confirms that the edge algebras on the path are the endomorphism algebras of the complexes the Hom test examined (so E-159's "End(T) matched to the child by Hom dimensions" caveat is discharged for these 13 edges), not that the complexes generate K^b(proj). Isomorphism is label-preserving; no vertex permutation was tried. Refuted by: a path edge where the Groebner ideal contains 1 or the arrow counts differ.

## Evidence
- 13 edges, `toolsmith_endt_out.txt` (section path13): every line "dims same, arrows same, ... iso". Max dim Hom(T_i, T_j) = 2 on 6 edges (needs the correction terms).
- Power of the decision procedure (section perturb): for each non-monomial relation of c on the 13 edges (22 relations) replace the binomial by its first term (a monomial End(T) does not satisfy): 22 of 22 rejected ("NO", 1 in the ideal). Doubling one coefficient is rejected in 10 of 22 (the rest absorbed by arrow scalars, as expected when the torus system has no cycle through that relation; the 3 edges whose torus system has a cycle, lna9 R1, F2, F7, reject 10 of 10).
- Selftest (section selftest): all 58 gate-admitted vertex-preserving steps on the first 12 LNAs of length 5 and duals: dims and arrows same, relations vanish or the scalar system is consistent; all J = 0, tiltingPlus true.
- Controls missing: no step with a different (wrong) child of equal dims and arrows was available, so rejection power is shown by perturbed relations only, not by a wrong algebra from the walk.

## Reproduction
```
DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 1 20000 /tmp/tsm/c1.pkl 100   # 513 s; rebuilds the pickle (nothing survives rounds)
.venv/bin/python workshop/rounds/053/toolsmith_endt_run.py path13     # ~2 s each; also: selftest | tilt | perturb | fail
```
Library part `workshop/rounds/053/toolsmith_endt.py` (class `TiltEnd`); recorded output `workshop/rounds/053/toolsmith_endt_out.txt` (11 KB).

## Prior record
No quiver-level End(T) result in `research/` (grep "End(T)", "quiver level"): E-159 states the caveat "End(T) is matched to the child by Hom dimensions, not as a quiver with relations"; this discharges it for 13 edges. Nothing in RETRACTIONS.md bears on it.

## Code changed
None in the library. New: `toolsmith_endt.py`, `toolsmith_endt_run.py`, `toolsmith_endt_out.txt` (round 053). No tests needed (no library file touched); self-checks above.

## Next
- skeptic: replay with an independent End(T) (e.g. via the Gabriel quiver of modules over a different convention) is of limited value; more useful is a test of generation (T generates K^b(proj A)): at an irreducible mutation it holds iff P_v -> add(T/P_v) is a left approximation with cokernel as stated (AI 2.31/2.32 hypothesis) -- the theorist's item.
- toolsmith: extend to class 2, the 3 E-158 paths and the 25 E-155 paths (only `/tmp/tsm` pickles needed, ~10 min each); parallel-arrow steps need a matrix-valued arrow identification (8 of 16 failing steps now undecided).
- Because End(T) agrees at J != 0 steps, the failing steps differ from the 13 only by Hom(T, T[-1]) != 0 (E-126/E-159); the discriminating object is the extension group, not the endomorphism algebra.
