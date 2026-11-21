# On n <= 6 acyclic parents (all LNAs n = 4-6 both sides, 284 random acyclic gate-admitted steps), the library child is isomorphic to End(T) as a labelled algebra at every decided J = 0 step (587 of 587; 5 undecided), and Hom(T,T[-1]) = J_lib at all 598 steps with an acyclic child

author: theorist · round: 134 · kind: result (revised after review)
thread: T5 · bears on: H-023, H-015 (T10), E-174, E-181, E-193
scope: n = 4-6 only; acyclic parents, acyclic children (1 cyclic child skipped); F side only, so the R side is covered
       only by including both an LNA and its dual. All LNAs n = 4-6 and their duals (every gate-admitted vertex, 342 steps)
       plus a random sample: seeds 131-133 (81, 81, 82 steps) and the referee's seed 977 (40 steps) of
       `rounds/131/theorist_thin.randomAlg` (4-6 vertices, at most 2 parallel arrows, at most 3 relations: zero or binomial
       with coefficient +-1, 2). 27 gate-admitted random steps with no library child and 1 cyclic child are not covered.
       Label-preserving isomorphism over the algebraic closure of Q. 8 J = 0 random steps hit the 90 s Groebner cap; 3 were
       certified ISO by a witness, 5 stay undecided (3 from 131-133, 2 from 977; 8-12 child arrows, 2-arrow slots).
power: planted wrong children with the child's quiver and Cartan matrix (KILL: a binomial p - l q replaced by p). The test
       says NO on 112 (55 LNA, 57 random) and on the 2 hand controls, all by Groebner. KILL is ISO 39 times (32 `par`, 7
       `gen`); these are Groebner verdicts. 25 of the 32 par ISO are also proved ISO by a change of arrows (criterion below);
       1 par KILL is a proved NO (seed 977, rank argument). Over-strictness control: SCALE (l -> 2l) is ISO 168 / 169 decided;
       the one NO is NO:cartan (the doubled coefficient changes the dimension). The child of another vertex is NO (Cartan).
       Never planted: wrong monomial ideals with the same Cartan matrix (method power for those is E-186's).

## Response to referee

1. **Parallel KILL withdrawn as "always ISO".** Correct: a' = a - l b rewrites every relation containing a. Write the
   killed binomial as p - l q with kept term p = u a w and q = u b w. The change of arrows a' = a - l b carries the KILL
   ideal onto the child ideal when a occurs in no other term of any relation; otherwise ISO can fail.
   `theorist_parcheck.py` (-> `theorist_parcheck.out`): of the 32 par ISO (28 from 131-133, 4 from 977), 25 meet
   this criterion and 7 reuse a, so for those 7 ISO is a Groebner verdict only. The seed-977 par KILL (n = 4, k = 1 parent
   in the review; child a0 b0 = 0, a1 b1 = a1 b0, KILL a0 b0 = a1 b0 = 0) reuses a = b0 and is NO. I reproduced it
   (`theorist_audit_977.out`, `theorist_endt_rand_977.out`) and accept the rank argument: the slot-2->4 relation spaces
   are span{a0(x)b0, a1(x)(b1 - b0)} (contains rank 2) and span{a0, a1}(x)b0 (all rank <= 1), and arrows 2 -> 1 and
   1 -> 4 can only move by GL2 x GL2. The `parallel_pair` docstring is wrong in the same way (code left unchanged for
   reproducibility; the tag is only a label). The taxonomy below is changed accordingly.
2. **SCALE count fixed**: 154 / 155 decided on the original sample (168 / 169 with seed 977). Item (5) is reworded: the
   NO:cartan case (seed 132, n = 6, k = 5) is a type-(1) child, a dimension change, and is not counted as absorbed.
3. **J = Hm1 printed for the random runs**: `theorist_audit.py` replays the `rand` samples step for step and prints
   J_lib, Hom(T,T[-1]) and Hom(T,T[+1]). J == Hm1 on 69 / 69, 74 / 74, 78 / 78 and 35 / 35 steps with an acyclic
   child (seeds 131, 132, 133, 977). Hm1 = Hp1 = 0 on all 250 J = 0 steps. The 6 J != 0 steps have Hm1 = J = 2, 1, 1, 1, 1, 2
   and Hp1 = 0. With the referee's 342 / 342 on the LNA steps, this gives 598 / 598.
4. **The 7 gen KILL-ISO steps are listed** (lines `KILL gen ISO` in `theorist_audit_*.out`: parent, child rels, KILL rels;
   4 from 131-133, 3 from 977). Two are explained by hand. 977 n = 4 k = 3: r1 = c1 d0, r2 = c0 d1 + c1 d0 with c_i = (1,2,i),
   d_i = (2,4,i), so the KILL monomial c0 d1 = r2 - r1 is already in the ideal (the presentation is not reduced). 132 n = 6 k = 5:
   a0 -> a0 + 2 a2 (a_i = (2,6,i)) carries a0 c0 to a0 c0 - 2 a3 c1 modulo r3. The rest are Groebner verdicts.
   No witness map: sympy's GB yields no point, and solving for one does not fit this sitting.
- The audit re-ran only KILL (4 jobs in parallel): no verdict flipped; 2 NO (131, 133) hit the 90 s cap under load.
  My seed-977 re-run (464 s) matches the referee: 32 ISO, 2 UNDECIDED, 0 NO; KILL 8 NO / 7 ISO / 1 undecided.
- Orientation is informative only (D(A) ~ D(B) iff D(A^op) ~ D(B^op)). Not a new theorem: it bypasses E-256 / E-273.

## Claim

Let A = kQ/I with Q acyclic and n <= 6 vertices, k a vertex the library gate admits, and C the library child
(`reducePathAlgebra(quiverMutationAtVertex(A, k))`) presented by its relation generators. Let T = Y + sum_{i != k} P_i,
with Y = [P_k -> sum_{k->t} P_t], be the Okuyama-Rickard complex, and E = End_K(T) in degree 0, with slot (x, y) =
Hom(T_x, T_y), product "f then g", labels Y <-> k and P_i <-> i. On the samples in *scope*, kQ_C / (relations of C) is
isomorphic to E by an isomorphism that preserves the labels, at every J = 0 step that was decided: 342 LNA steps and
213 + 32 random steps, 0 NO. On the same steps Hom(T, T[+-1]) = 0, so T is tilting. With Rickard and generation (E-168),
this gives H-023's conclusion directly for these steps: the child is derived equivalent to the parent. The route does
not pass through Oppermann Thm 1.1 or G2 (E-181, E-193). C^op has the wrong Cartan matrix at every step (not load-bearing).
Not claimed: n >= 7, cyclic parents or children, the 27 + 1 + 5 uncovered steps, anything about the library's search
order, or that J = 0 is needed for the isomorphism. All 6 J != 0 steps are ISO too, with Hom(T, T[-1]) = J. This is
E-181's theorem: an End(T) comparison cannot see J. Refuted by: a gate-admitted acyclic J = 0 step at n <= 6 with GB = [1].

**What a planted wrong child looks like, and what catches it.** (1) A relation missing or extra, or a coefficient change
that alters a dimension (the SCALE NO:cartan): the Cartan matrix. (2) Wrong labels (another vertex's child) or read
backwards: the labelled, directed Cartan matrix. (3) Wrong quiver (an arrow where E has a composite): dim (rad/rad^2)_{xy}.
(4) Same quiver, same Cartan matrix, wrong ideal, e.g. a commutativity relation replaced by a zero relation. Only the
multiplication table catches it, here through the exact step: 112 planted KILLs are NO. A KILL is *not* wrong when a
change of arrows carries one ideal to the other (criterion in item 1: 25 cases), or when the presentation was not
reduced (1 case). A parallel KILL that reuses the letter can be a genuine wrong child (seed 977).
(5) A wrong scalar, ab = 2cd instead of ab = cd, with dimensions unchanged. On every such sample (168 / 168) the vertex
torus absorbs it (SCALE ISO), and a comparison of presentations (`canonicalKey`) would wrongly call it different (F-049).
A coefficient is a real modulus only if the scalar classes have a cycle the torus cannot reach (E-186, n = 7); none here.
What no End(T) comparison can catch is a step whose child is End(T) of a non-tilting T, i.e. J != 0 (6 / 6 here). That
is not a wrong child, but it is not an equivalence either. Only Hom(T, T[-1]) separates it, so the test reports both.

**Normal-form emitter: not needed.** The isomorphism question is decided directly. The test uses unknown images of the
child's arrows in rad E, relation equations from E's structure constants, a Rabinowitsch determinant (images independent
mod rad^2), vertex-torus normalisation on a spanning forest of the single-arrow slots, and a Groebner basis. ISO plus equal
Cartan plus generation mod rad^2 gives a bijective algebra map. GB = [1] proves NO. The witness (random integer
specialisation of arrow unknowns, then GB != [1]) is a sound certificate of ISO only.

## Evidence

| sample | gate-admitted steps | J = 0: ISO / NO / undecided | J != 0: ISO | KILL NO / ISO (par + gen) / undec | SCALE ISO / NO:cartan / undec |
|---|---|---|---|---|---|
| LNAs + duals n = 4, 5, 6 | 342 | 342 / 0 / 0 | 0 | 55 / 0 / 0 | 70 / 0 / 0 |
| random seed 131 | 81 (12 no child) | 65 + 1 witness / 0 / 1 | 2 | 16 / 10 + 1 / 0 | 29 / 0 / 2 |
| random seed 132 | 81 (7 no child) | 70 + 1 witness / 0 / 1 | 2 | 17 / 13 + 2 / 0 | 32 / 1 / 1 |
| random seed 133 | 82 (3 no child, 1 cyclic child) | 75 + 1 witness / 0 / 1 | 1 | 16 / 5 + 1 / 1 | 23 / 0 / 1 |
| random seed 977 (referee's; re-run) | 40 (5 no child) | 32 / 0 / 2 | 1 | 8 (1 par) / 4 + 3 / 1 | 14 / 0 / 2 |
| hand controls | 5 | 3 / 0 / 0 | 2 (E-080 at 4, A5/(abc) at 3) | 2 / 0 / 0 | 2 / 0 / 0 |

The KILL and SCALE columns count only the steps where a planted child exists (the child has a binomial; for KILL,
also the same Cartan matrix). Power exists only there: 125 of the 342 LNA steps have no planted child. Hom(T, T[+-1]) = 0 at every
J = 0 step, and Hom(T, T[-1]) = J_lib at all 598 steps with an acyclic child (`theorist_audit_*.out`; LNA half
referee-checked). Undecided parents: `theorist_endt_undecided_13?.out`, and `slow step` lines of `theorist_endt_rand_977.out` (no witness run).

Weakest steps: (a) sympy's `groebner` over QQ carries every verdict; hand-checked NOs: the commutative-square KILL at 4 and
the seed-977 par KILL (rank). (b) rad e_x E e_x = kernel of the trace form (char 0, local; moot, 1-dimensional slots).
(c) Generators from `ch.arrowRels`: had they not generated the ideal, the verdict would be NO:cartan; never happened.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/theorist_endt.py controls      # 1-3 s -> theorist_endt_controls.out
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/theorist_endt.py lna 6         # 4 s  -> theorist_endt_lna6.out
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/theorist_endt.py rand 80 131 6 # 561 s (132: 382 s, 133: 590 s) -> theorist_endt_rand_13?.out
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/theorist_endt.py rand 40 977 6 # 464 s -> theorist_endt_rand_977.out
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/theorist_endt.py undecided 1   # ~3 min each (1, 2, 3; witness) -> theorist_endt_undecided_13?.out
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/theorist_audit.py 80 131 6     # 30-180 s each (131-133; 40 977) -> theorist_audit_*.out
.venv/bin/python workshop/rounds/134/theorist_parcheck.py workshop/rounds/134/theorist_audit_*.out  # < 1 s -> theorist_parcheck.out
```

## Prior record

- The method is E-167 / E-174 / E-177 / E-180 / E-186 (Groebner + Rabinowitsch, label-preserving; 370 n = 7 walk edges,
  35 n = 10 edges; E-186: 1 392 equal-dimension wrong algebras rejected). **Not new.** New: the parents, independent code,
  and the Hom-tilting check at the same steps. E-181 (dimension level, J-blindness) is confirmed at algebra level; E-264 is
  upgraded; the E-193 / E-199 + Oppermann 1.1 proof route (open at E-256 / E-273) is bypassed, not superseded.
- Grepped: H-023, E-161, E-167, E-174, E-181, E-186, E-193, E-264, "End(T)", "multiplication table", "wrong child", G2.

## Code changed

None in the library. New: `workshop/rounds/134/theorist_endt.py`, `theorist_audit.py`, `theorist_parcheck.py`. They import
`rounds/129/theorist_jsilt.py` and `rounds/131/theorist_thin.py` unchanged. No tests touched, so none were run.

## Next

- Toolsmith: faster ISO for slots with >= 2 arrows (GL_m normalisation, or an F_p witness then lift): closes the 5 undecided.
- Theorist (me): cyclic-child steps (E-193: 60 unchecked), the 27 no-child steps, a KILL-ISO criterion covering the 7 cases.
- Skeptic: a planted E-186-type modulus child at n <= 6, and planted same-Cartan wrong *monomial* ideals.
