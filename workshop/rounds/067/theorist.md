# g = gd (proof, any quiver, true paths; not computed by any code); step-7 kernel = K when D(i) is acyclic, including cycles upstream of D(i) (7 389 pairs, n = 7..9 total vertices, 0 differences; checks the repo code, not the mathematics)

## Response to referee

1. thruI counter: fixed in `theorist_cycup.py`; it now counts pairs (Q, i) where some relation starting on the attached cycle visits i. Re-run: 7 389 pairs, thruI = 1 692, all other numbers unchanged (K != 0 727, g = gd = 859, 0 differences, `ALL EQUAL, g = gd`). So the "through i" claim is checked in 1 692 pairs.
2. 6 080 corrected to 5 840 (the referee's re-run is right; my text was a typo; re-run of cycoff gives on = 5840, off = 608, offk = 1).
3. Title, claim and scope narrowed as proposed: g = gd is a proof only (no code computes true-path Z_1); the upstream-cycle test is vacuous in principle by (C) and checks the repo code, not the mathematics; "n = 5..9" replaced by 5..7 old vertices plus 2 cycle vertices.

author: theorist · round: 067 · kind: proof
thread: T5 · bears on: H-015 (premise J = 0), E-193, E-189
scope: proof for arbitrary quivers with TRUE paths (no acyclicity, parallel arrows allowed, relations of length >= 2, any field); the repo-side identification only for parents with D(i) acyclic; tests: 5..7 old vertices plus an attached 2-cycle (7..9 vertices in total); the upstream-cycle test is vacuous in principle (by (C) the cycle is invisible to step 7), so it checks the repo implementation only, and the true-path statement (A) is untested (proof only), random multigraphs, rational coefficients, 12 seeds. NOT covered: a cycle meeting D(i) (downstream of i or through i).

## Claim

Notation (E-193): Z_1(i,v) = cycles in the degree-1 part of Oppermann's dg algebra from i to v, spanned by words (q, r, p) (q: i -> s(r), p: t(r) -> v) with d(q r p) = q rho_r p; B_1 = boundaries of degree-2 words; red(z) = the part of z with q empty, as a vector in {(r, p)}; K(v) = red(Z_1(i,v)); g(v) = dim Z_1/(B_1 + Z_1 rad + rad Z_1); gd(v) = dim K(v)/(sum_a K(w) a + red B_1(v)).

(A) g(v) = gd(v) for every v, for every quiver and every relation set, if "paths" are true paths. So each minimal generator x has red(z_x) independent modulo decomposables.
(B) red(B_1 e_i) = rbar . I (all (r, s) with r a relation at i and s in I(t(r), v)), and K(v) = { z_0 : (d z_0)/a in I for every arrow a out of i }, again for true paths.
(C) The repo's step 7 computes (B) correctly whenever the full subquiver D(i) of vertices reachable from i (including i) is acyclic: paths and ideal membership below i only see D(i). An oriented cycle upstream of or disjoint from D(i) never enters step 7.

It does NOT claim: anything when a cycle meets D(i) (there the repo's vertex-simple `allPathsBetween` / `isInIdeal` are wrong, E-193, and true Z_1 may be infinite); anything about H_1(M) beyond the generators (see "J = 0" below); the cancellation step of Thm 1.1 or the left/right convention. Refuted by: one pair (Q, I, i) with true-path g != gd, or D(i) acyclic and repo kernel != K.

## Evidence

Proof of (A). red is linear Z_1 -> K and surjective by definition of K. It maps B_1 onto red B_1, Z_1 rad onto K rad (red(z a) = red(z) a), and rad Z_1 to 0 (a z has q nonempty). So red induces a surjection Z_1/(B_1 + Z_1 rad + rad Z_1) -> K/(K rad + red B_1), g >= gd. Injective: let red z lie in K rad + red B_1; choose y in Z_1 rad + B_1 with red y = red z. Then z' = z - y is a cycle with red z' = 0, i.e. every word of z' has q nonempty. Write z' = sum_a a z'_a over the first arrow a of q (words of a z'_a are distinct for distinct arrows a, parallel arrows included). d z' = sum_a a . d(z'_a) = 0 in the free path algebra, and the summands lie in different a-blocks, so d(z'_a) = 0 for each a: z'_a in Z_1(t(a), v), z' in rad Z_1. Hence z in B_1 + Z_1 rad + rad Z_1. QED. Weakest step: none; it uses only that left multiplication by a distinct arrow is injective on the free path algebra and that words are a basis.

Proof of (B). B_1 is spanned by d(q r2 m r1 p) = (q rho2 m) r1 p - q r2 (m rho1 p) (signs as in theorist_g2exact.py). In the first term the tail q rho2 m is non-empty because relations have length >= 2, so red kills it. In the second, the tail q is empty iff s(r2) = i, and the term is (r2, m rho1 p): as (m, r1, p) run over all true paths, m rho1 p spans I(t(r2), v). Cycle condition: z = z_0 + sum_a a z'_a is a cycle iff d z_0 = - sum_a a d(z'_a); d z_0 has every term beginning with an arrow, d(z'_a) ranges over all of I(t(a), v), so z_0 in K iff (d z_0)/a in I for all a. QED. (Lemma L of E-193 is the part B_1 -> rbar . I; E-193 stated it "valid where no path from t(r) returns to i": that was the vertex-simple restriction, not part of the mathematics.)

Proof of (C). Every path or ideal element occurring in (B) goes from a vertex of D(i) to a vertex of D(i) (between x and y both below i), so only D(i) matters; if D(i) is acyclic, true paths there are vertex-simple and the repo's objects equal the true ones; E-193's acyclic-parent equality then applies verbatim (it used only these paths).

Test (cycle off D(i)). `theorist_cycup.py`: E-193's random acyclic instances (with relations) plus an oriented 2-cycle c1 <-> c2 attached UPSTREAM (arrows c -> old vertices, extra relations starting on the cycle and running through old vertices, so through i for some i). For every old i with an out-arrow, the repo's raw step-7 kernel vs red(Z_1 e_i) as subspaces, and g vs gd:

| pairs | with K != 0 | multi-term K vector | parallel arrows in Q | differences K | sum g | sum gd | g != gd | exceptions | #relations from i* after `_minimalGenerators` = g |
|---|---|---|---|---|---|---|---|---|---|
| 7 389 | 727 | 411 | 6 021 | 0 | 859 | 859 | 0 | 0 | 7 389 of 7 389 |

(In 1 692 of the 7 389 pairs a relation starting on the cycle passes through i, counted by `thruI`, now incremented. The count equals g in all 7 389, including pairs whose child is cyclic; E-193 had 7 478 of 7 479.) `theorist_cycoff.py` (random back edges, classifies 'on'/'off'): 5 840 pairs 'on' (a cycle meets D(i)) are not evaluated (the vertex-simple Opp is meaningless there), 608 pairs 'off' agree, but only 1 has K != 0, so that run has almost no power; the powered test is cycup.

## What this means for J = 0

1. G2 (E-193) is now a theorem for every parent with D(i) acyclic, with the step-7 generator COUNT proved to be g (the count of degree-2 arrows i -> v of Q^(2)) up to the redundant Lemma-L elements `_minimalGenerators` drops (that dropping is observed, not proved). The acyclic-parent restriction in E-193/STATE shrinks from "no oriented cycle anywhere in the parent" to "no oriented cycle meeting the descendants of the mutated vertex i". The toolsmith scan of printed paths needs only this condition, not global acyclicity; since the repo's mutation preserves everything above i, many cyclic parents (cycle upstream) are already inside the proved range.
2. g = gd closes one candidate source of J != 0: no degree-2 generator of kM e_i has zero image in H_0 (gd = g), i.e. there is no "x with red(z_x) = 0". It does NOT prove J = 0: J lives in H_1(M), whose classes need not be generators; (A) is a statement about Z_1/(decomposables) and says nothing about relations among generators. The premise still rests on the per-step linear-algebra check J_i = 0 (E-099/E-171), cancellation (E-189) and the left/right convention.
3. For cyclic parents with a cycle in D(i) no result here; the true-path version of (A), (B) is proved but computing it needs a repo that enumerates true paths modulo I (finite only when the algebra is finite-dimensional).

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/067/theorist_cycup.py 12 250   # 48 s, last line ALL EQUAL, g = gd
timeout 10m .venv/bin/python -W ignore workshop/rounds/067/theorist_cycoff.py 12 300   # 3 s, low power (1 pair with K != 0)
```
(`theorist_cycup.py` imports `workshop/rounds/065/theorist_g2exact.py`.)

## Prior record

E-193 proved K = Ker for acyclic parents and recorded "g = gd observed, unproved" and Lemma L "valid where no path from t(r) returns to i"; E-189 (generators in degree <= 2). New: (A) the proof of g = gd, (B) Lemma L and the K description without acyclicity, (C) the D(i) criterion, and the upstream-cycle test. Nothing in `research/RETRACTIONS.md` concerns this (grepped g = gd, Lemma L). The proofs are short enough that a referee should hand-check (A) first.

## Code changed

None in `quivermutation/` (the CyclicQuiverWarning already present was suppressed with `-W ignore` in the runs). Two new scripts in `workshop/rounds/067/`. No test files touched.

## Next

- toolsmith: rescan printed paths (E-157/E-160/E-163/E-169) with the weaker condition "a cycle meets D(i) of the mutated vertex" at each step, instead of "parent has a cycle".
- skeptic: hand-check (A) on the smallest g = 1 pair with a parallel arrow, and the sign of B_1's second term.
- theorist next: a cycle in D(i) with a finite-dimensional algebra (cyclic Nakayama i -> cycle): compute K via true paths mod I and compare with the repo; then H_1(M) beyond generators (is J = 0 equivalent to a per-vertex statement about Z_1 e_i / (B_1 + decomposables) in higher degree?).
