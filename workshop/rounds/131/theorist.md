# At a gate-admitted vertex v, every i != v with dim e_i A e_v <= 1 has J_i = 0 (E-128 L1, any field, any relations, parallel arrows allowed); with all such dims <= 1 and I admissible (R nilpotent mod I, a separate hypothesis), even on a cyclic quiver every nonzero path into v is vertex-simple, so the library's simple-path gate is the exact gate there; J != 0 at an admitted step needs dim e_i A e_v >= 2 and a non-monomial piece of I (checked per row: 29 of 29 and 92 of 92 rows)

author: theorist · round: 131 · kind: proof (core step is a rediscovery of E-128 L1)
thread: T10 / T5 · bears on: H-032, H-023, E-337, E-333, E-130, E-171, E-128, E-131, E-134, E-208, E-310, E-329
scope: proof for any finite-dimensional A = kQ/I (k any field); part (c) also needs the arrow ideal R nilpotent mod I (I admissible), which does not follow from finite dimension plus I in R^2. Checks: 600 + 600 random acyclic algebras on 4-7 vertices (zero relations and two-term commutativity relations with coefficients ±1, 2, up to 2 parallel arrows), 1 600 admissible vertices each; (d) checked per row on the 29 gate-admitted J != 0 rows of that sample and on the 92 of a 2 000-algebra sample (seed 7); 600 random cyclic monomial algebras on 3-5 vertices (1 739 thin vertices); E-329's bag instance as a record check of (c); acyclicity of the step quivers in 13 F-041 balls (n = 10: 11 pairs, n = 11: 2 pairs, 6 070 admitted steps). Not covered: non-monomial cyclic algebras (lemma (c) is proved but not computed on them); whether the library's `isInIdeal` is exact on cyclic quivers (E-329 says it is not); anything about child = End(T) (H-023).
power: n/a for the proof. Computational check: a planted failure for each dropped hypothesis is caught. (1) Drop d <= 1: E-080 square at v = 3 (d = 2), gate True, J = 1, library `isTiltingAtVertex` False. E-318's parent at v = 3 (d = 2) is the same: gate True, tilting False, while its six other admitted vertices have d <= 1 and are tilting. (2) Drop the gate: 1 -> 2 -> 3 with ab = 0 at v = 2 (d = 1), J = 1. (3) Drop thinness in the cyclic step: 380 thick vertices in the cyclic sample have a nonzero path into v that is not vertex-simple.

## Response to referee

1. **(c) weakest-step remark: corrected.** The referee is right: finite dimension plus I in R^2 does not give R^m in I (k[x]/(x^2 - x^3) has dim 3 and x is not nilpotent, since the factor k[x]/(1 - x) sees x = 1). The proof of (c) uses only C^m in I, so nilpotence of R mod I (I admissible) is now stated as a separate hypothesis in the title, Scope and Claim (c), and the false sentence is replaced. The result of (c) is unchanged.
2. **(d) per row: computed.** New `theorist_dcheck.py rows` replays the sweep (same RNG stream; the row counts reproduce 29 and 92) and tests, at each gate-admitted J != 0 vertex, whether some i with d_i >= 2 and some out-arrow beta have e_i I e_t(beta) not spanned by paths. Result: 29 of 29 (`rows 600 131`) and 92 of 92 (`rows 2000 7`, the referee's wider sample), also when i is restricted to the i with J_i != 0. The test is not constant (it fires at 95 of 287 and 289 of 969 admitted J = 0, d >= 2 rows). Title and Scope now say (d) is checked per row.
3. **E-329 bag: added** as a record check of (c) (`theorist_dcheck.py bag`, independent engine): d_5 = 2 (g, bag), so v = 1 is thick, and the only refusing path, bag, is non-simple. Exact gate False, matching E-329.

## Claim

Setting: A = kQ/I is finite-dimensional, I lies in the arrow ideal R, and v has no loop and at least one arrow out. For i != v put d_i = dim e_i A e_v. J_i = {x in e_i A e_v : x beta = 0 in A for every arrow beta out of v}, counted by arrow, so parallel arrows count twice (Ladkani 2.3(c); E-128; `isTiltingAtVertex`). The gate G(v) (`procedure.isMutable`) says that no path p: i ~> v with i != v is nonzero in A while p beta = 0 for every out-arrow beta. The library enumerates only vertex-simple p; call that version G_s(v).

(a) G(v) implies that J_i contains no nonzero path image, so dim J_i <= d_i - 1 whenever d_i >= 1. In particular **d_i <= 1 together with G(v) gives J_i = 0**, one i at a time. This is E-128 (L1), proved in round 034. It is not new.
(b) Conversely, J = 0 implies G(v). So at a thin vertex (all d_i <= 1) G(v), J = 0 and tilting are the same condition. Tilting is equivalent to J = 0 by E-130 / E-171 (AI 2.31, 2.32(b), any finite-dimensional A, v loopless).
(c) (New.) Assume R is nilpotent mod I (I admissible; this is a hypothesis in its own right, not a consequence of dim A < infinity) and every d_w <= 1 (w != v). Then every path from i != v to v that is nonzero in A is vertex-simple. So at a thin vertex G_s(v) = G(v) also on a cyclic quiver, and (a) holds for the gate as the library states it.
(d) Shape of every exception. If v is gate-admitted and J != 0, then for some i != v and some out-arrow beta both of the following hold: d_i >= 2, and e_i I e_t(beta) is not spanned by paths. A generator of J_i is x = sum c_p p with at least two nonzero c_p in every path basis of e_i A e_v, x is not in I, and x beta is in I for all beta. The minimal case is d_i = 2 and x = p - lambda q with (p - lambda q) beta in I: the E-080 / E-318 shape. Hence for a monomial I, G(v) and J = 0 are equivalent at every d (already in the record: E-208, E-310).

Not claimed: (1) that a thin step is a derived equivalence. Thinness gives T tilting, and child = End(T) is H-023, which is open. (2) Anything about the library's ideal membership on cyclic quivers (E-329). (3) That the F-041 balls must be thin. That is E-337's measurement, not a theorem.
Refuted by: an algebra kQ/I and a vertex v with G(v) true, some d_i <= 1 and J_i != 0. For (c), a nonzero non-simple path into v when all d_w <= 1 (w != v).

## Evidence

**Proof of (a), (b).** Paths from i to v span e_i A e_v. If d_i >= 1, a nonzero element of e_i A e_v is a combination of path images, so some path p has p != 0 in A. If d_i = 1 then e_i A e_v = k·p. G(v) gives an out-arrow beta with p beta != 0. So x -> (x beta)_beta is nonzero on a one-dimensional space, hence injective, and J_i = 0. For general d_i, G(v) says that J_i meets the set of path images only in 0, and the path images span. So J_i is a proper subspace, dim J_i <= d_i - 1. For (b): a path that violates G(v) is a nonzero element of J_i.
Hypotheses used: the paths span e_i A e_v (any I inside R); the gate's quantifier ranges over *all* nonzero paths into v; out-arrows are counted as arrows; v has an out-arrow (with none, J_i = e_i A e_v). Not used: acyclicity, admissibility, monomiality, the field, the out-degree, the key guard, or that the parent lies on a walk from an LNA.

**Proof of (c).** *Weakest step: the nilpotence argument. It needs C^m in I for some m, which holds when R is nilpotent mod I (I admissible). That is an extra hypothesis: finite dimension plus I in R^2 does not give it (one vertex, one loop x, I = (x^2 - x^3) in R^2: A = k[x]/(x^2(1 - x)) has dim 3 but x^m != 0 for all m). (An earlier version of this sentence said the opposite; it was false, see Response item 1.)* Let p: i ~> v be nonzero in A and not vertex-simple. Let w be the first vertex that p visits twice, and write p = p1 C p2 with C a nontrivial cycle at w.
- If w != v: then p2 and C p2 are both nonzero in e_w A e_v, because each is a factor of p. If C p2 = lambda p2, then C^m p2 = lambda^m p2. C^m = 0 for large m, so lambda = 0, but then C p2 = 0, a contradiction. So the two are independent and d_w >= 2. This covers w = i.
- If w = v: write p = p1 Z with p1: i ~> v the first arrival at v and Z a nontrivial cycle at v. Then p1 and p1 Z are both nonzero, and by the same argument they are independent, so d_i >= 2.
Either case contradicts thinness. Then G_s(v) = G(v), and (a) applies.

**Proof of (d).** By (a), J != 0 forces some d_i >= 2, and no nonzero element of J_i is a path image, so its support has size >= 2 in every path basis. Suppose e_i I e_t(beta) were spanned by paths for every beta. Then path images (nonzero paths) form a basis of e_i A e_t(beta), distinct p give distinct p beta, and x beta in I would force p beta in I for each p in the support of x. The gate would then refuse each such p, a contradiction. (The statement in the Claim is for I monomial in the degrees (i, t(beta)), which is all the argument uses.)

**Per-row check of (d)** (`theorist_dcheck.py rows`, same RNG stream and engine as the sweep). At each gate-admitted J != 0 vertex: (d1) some i with d_i >= 2 and some out-arrow beta with e_i I e_t(beta) not spanned by paths (fully reduced echelon basis not all single paths); (d2) the same at an i with J_i != 0. `rows 600 131`: 29 rows, (d1) 29, (d2) 29. `rows 2000 7`: 92 rows, (d1) 92, (d2) 92. The test is not constant: at admitted J = 0 rows with d >= 2 a non-monomial piece is present at 95 of 287 and 289 of 969 (informative; (d) says nothing about these rows).

**Record check of (c): E-329's bag** (`theorist_dcheck.py bag`; a:4>5, b:5>4, g:5>1, x1:1>2, x2:2>3, relations ab, bagx1, v = 1). Independent engine: d_5 = 2 (basis g, bag), d_4 = 1, d_2 = d_3 = 0, so v is thick. The only refusing nonzero path is bag, and it is non-simple (visits 5 twice); the exact gate is False. This is the w != v case of the proof of (c) on a recorded instance, and it is why the simple-path library gate says True there.

**Computation** (`theorist_thin.py`, an independent exact engine: it builds e_i I e_j from u r w and row-reduces over Q; library maths is used only for the comparison):

| sample | (min(dmax,2), gate, J != 0) : count | library disagreements (gate or tilting) | dim J_i > d_i - 1 at admitted v |
|---|---|---|---|
| 600 acyclic, mixed relations | (0,T,F) 1060 · (1,T,F) 476 · (1,F,T) 169 · (2,T,F) 287 · (2,T,T) **29** · (2,F,T) 89 | 0 | 0 |
| 600 acyclic, monomial | (0,T,F) 1100 · (1,T,F) 515 · (1,F,T) 166 · (2,T,F) 250 · (2,T,T) **0** · (2,F,T) 94 | 0 | 0 |

No row is (<=1, T, T), and no row has gate F with J = 0. Gate-admitted J != 0 occurs only at d >= 2 with a non-monomial I, as (a) to (d) say. Cyclic monomial sample (part (c)): 1 739 thin vertices, path-simplicity held at all 1 739. At thick vertices, 380 have a nonzero non-simple path into v (the planted failure of (c)).

Planted / record checks (`theorist_thin.py planted`): E-080 square at v = 3: gate True, dmax 2, J 1, library (True, False). Gate dropped (ab = 0, v = 2): gate False, d 1, J 1. Monomial square: gate False, J 1. E-318 parent: v = 3 has d = 2 and is not tilting; v = 1, 2, 4, 6, 8, 9 have d <= 1 and are tilting.

F-041 balls, acyclicity (`theorist_thincyc.py`, reusing round-130 `toolsmith_f041calls.py`): n = 10 pairs 0, 1, 2, 3, 50, 100, 150, 200, 250, 302, 303 and n = 11 pairs 5, 131, depth 3 with duals, option on. 6 070 admitted steps, **0 on a cyclic quiver**, dmax 1, 0 refused. These pairs are acyclic, so E-337's simple-path d is the true d there, and by (a) "isTiltingAtVertex never False" in E-337 is a consequence of "d <= 1 at every step". It is not a second observation. The other 291 + 1 511 balls were not checked for acyclicity. By (c), the conclusion would survive a cyclic step there as long as membership in I was decided correctly.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/131/theorist_thin.py planted          # 3 s
timeout 10m .venv/bin/python workshop/rounds/131/theorist_thin.py sweep 600 131    # 2 s
timeout 10m .venv/bin/python workshop/rounds/131/theorist_thin.py cyclic 600 131   # 0.2 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/131/theorist_thincyc.py 10 0,303                        # 12 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/131/theorist_thincyc.py 10 1,2,3,50,100,150,200,250,302 # 44 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/131/theorist_thincyc.py 11 5,131                        # 13 s
cd workshop/rounds/131 && timeout 10m ../../../.venv/bin/python theorist_dcheck.py rows 600 131   # 0.3 s
cd workshop/rounds/131 && timeout 10m ../../../.venv/bin/python theorist_dcheck.py rows 2000 7    # 0.5 s
cd workshop/rounds/131 && timeout 10m ../../../.venv/bin/python theorist_dcheck.py bag            # < 1 s
```

## Prior record

- (a) is E-128 (L1) (round 034, my own persona, refereed). Maverick's thinness argument and E-337's explanation are a rediscovery of it: "J_i != 0 forces d_i >= 2" is printed in E-128's title. E-131 (walks: every J != 0 row has d = 2, dim J = 1) and E-134 (d = 3, J = 1 off-walk) are the data side.
- The monomial half of (d) (gate equals tilting for monomial I) is E-208 / E-310. The commutativity shape is E-080 / E-171's 5-vertex example.
- New: (c), the cyclic extension (thinness makes the simple-path gate exact). This answers, for thin vertices, E-130's caveat that "the repo gate enumerates simple paths only, so E-128's L1 is not guaranteed on cyclic quivers". Also new: the explicit hypothesis list, the converse (b), and the acyclicity count on 13 F-041 balls.
- For the T10 equivalence and the T5 premise: thin + gate gives J = 0, which gives T tilting (E-130, E-171, E-333). The step is then inside H-023's domain. Whether it is a derived equivalence still rests on H-023, not on this lemma.

## Code changed

None in the library. New scripts: `workshop/rounds/131/theorist_thin.py`, `workshop/rounds/131/theorist_thincyc.py`, `workshop/rounds/131/theorist_dcheck.py` (revision). No tests touched.

## Next

- Chair: record E-337's "because the balls are thin" as an instance of E-128 L1, with (c) as an addendum to E-130's cyclic caveat. The open part of T10's cost is only *where d >= 2 steps occur* (E-332 at depth 6, E-316's mixed walks), not why thin steps are tilting.
- Toolsmith: a cheap "d >= 2 and non-monomial in degree (i, t(beta))" pre-filter would let `refuseNonTilting` skip the full rank computation at all other steps (by (a) and (d), a valid shortcut).
- Theorist (next): H-023 at algebra level (End(T) multiplication table versus the child), which is the real gap behind both T5 and T10.
