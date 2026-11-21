# On n <= 6 acyclic parents (all LNAs n = 4-6 both sides, 244 random acyclic gate-admitted steps), the library child is isomorphic to End(T) as an algebra at every decided J = 0 step (555 of 555; 3 undecided), and the comparison rejects 104 of 104 planted same-Cartan wrong children it was given that are not isomorphic by a change of arrows

author: theorist · round: 134 · kind: result
thread: T5 · bears on: H-023, H-015 (T10), E-174, E-181, E-193
scope: n = 4-6 only; acyclic parents, acyclic children (1 cyclic child skipped); F side only, so the R side is covered
       only by including both an LNA and its dual. All LNAs n = 4-6 and their duals (every gate-admitted vertex, 342 steps)
       plus a random sample: 3 seeds x 81-82 gate-admitted steps of `rounds/131/theorist_thin.randomAlg` (4-6 vertices, at most 2
       parallel arrows, at most 3 relations: zero or binomial with coefficient +-1, 2). 22 gate-admitted random steps where
       the library returns no child (illegal relation, `child_of` None) and 1 cyclic child are not covered. Label-preserving
       isomorphism over the algebraic closure of Q. 6 J = 0 random steps hit the 90 s Groebner cap; 3 were then certified
       ISO by a witness, 3 stay undecided (8-13 child arrows, several 2-arrow slots).
power: planted wrong children with the child's quiver and Cartan matrix (KILL: a binomial p - l q replaced by p). The test
       says NO on 104 (55 LNA, 49 random) and on the 2 hand controls; e.g. the commutative square abd = acd at 4, where the
       KILL algebra is known to be non-isomorphic. KILL is also accepted as ISO 32 times: 28 when the two terms differ in one
       parallel letter (then a(b - l b') = 0 is a change of arrows, so ISO is correct) and 4 more where an isomorphism was
       certified. Over-strictness control: SCALE (l -> 2l) gives ISO 154 / 154 decided (NO:cartan 1, where the doubled
       coefficient changes the dimension). The child of a different vertex is NO (Cartan).

## Claim

Let A = kQ/I with Q acyclic and n <= 6 vertices, k a vertex the library gate admits, and C the library child
(`reducePathAlgebra(quiverMutationAtVertex(A, k))`) presented by its relation generators. Let T = Y + sum_{i != k} P_i,
with Y = [P_k -> sum_{k->t} P_t], be the Okuyama-Rickard complex, and E = End_K(T) in degree 0, with slot (x, y) =
Hom(T_x, T_y), product "f then g", labels Y <-> k and P_i <-> i. On the samples in *scope*, kQ_C / (relations of C) is
isomorphic to E by an isomorphism that preserves the labels, at every J = 0 step that was decided: 342 LNA steps and
210 + 3 random steps, 0 NO. On the same steps Hom(T, T[+-1]) = 0, so T is tilting. With Rickard and generation (E-168),
this gives H-023's conclusion directly for these steps: the child is derived equivalent to the parent. The route does
not pass through Oppermann Thm 1.1 or the G2 identification (E-181, E-193). The orientation is direct: C^op has the
wrong Cartan matrix at every step. Not claimed: n >= 7, cyclic parents or children, the 22 + 3 uncovered steps,
anything about the library's search order, or that J = 0 is needed for the isomorphism. All 5 J != 0 steps in the
sample are ISO too, with Hom(T, T[-1]) = 1, 1, 1, 1, 2. This is E-181's theorem: an End(T) comparison cannot see J.
Refuted by: a gate-admitted acyclic J = 0 step at n <= 6 whose Groebner basis is [1] (NO), e.g. on a re-run with other
seeds.

**What a planted wrong child looks like, and what catches it.** A wrong child can fail in five ways:
(1) A relation is missing or extra, so the dimension changes. The Cartan matrix catches this.
(2) The labels are wrong (the child of another vertex) or the child is read backwards. The labelled, directed Cartan
matrix catches this.
(3) The quiver is wrong (an arrow where E has a composite). dim (rad/rad^2)_{xy} catches this.
(4) Same quiver, same Cartan matrix, wrong ideal. For example, a commutativity relation is replaced by a zero relation:
this keeps every dimension and changes the rank of one multiplication map. Only the multiplication table catches it,
here through the exact step: 104 of 104 caught.
(5) A wrong scalar, ab = 2cd instead of ab = cd. This is not a wrong child: on every sample the vertex torus absorbs
it (SCALE ISO), and a comparison of presentations (`canonicalKey`) would wrongly call it different (F-049). A
coefficient could only be a real modulus if the scalar classes had a cycle the torus cannot reach. E-186 built such a
family at n = 7; none occurs here.
What no End(T) comparison can catch is a step whose child is End(T) of a non-tilting T, i.e. J != 0 (5 / 5 here). That
is not a wrong child, but it is not an equivalence either. Only Hom(T, T[-1]) separates it, so the test reports both.

**Normal-form emitter: not needed.** The isomorphism question here is decided directly; a normal form would have to
quotient presentations by Aut(kQ), which is exactly the SCALE / parallel-KILL freedom above. The smallest thing that does
the job is: unknown images of the child's arrows in rad E, relation equations from E's structure constants, a
Rabinowitsch determinant (images independent mod rad^2), vertex-torus normalisation on a spanning forest of the
single-arrow slots, and a Groebner basis. ISO plus equal Cartan plus generation mod rad^2 gives a bijective algebra map.
A Groebner basis equal to [1] is a proof of NO. The witness route (random integer specialisation, then a Groebner
basis not equal to [1]) is a sound certificate of ISO only.

## Evidence

| sample | gate-admitted steps | J = 0: ISO / NO / undecided | J != 0: ISO | KILL NO / ISO (par + certified) / undec | SCALE ISO / NO:cartan / undec |
|---|---|---|---|---|---|
| LNAs + duals n = 4, 5, 6 | 342 | 342 / 0 / 0 | 0 | 55 / 0 / 0 | 70 / 0 / 0 |
| random seed 131 | 81 (12 no child) | 65 + 1 witness / 0 / 1 | 2 | 16 / 11 / 0 | 29 / 0 / 2 |
| random seed 132 | 81 (7 no child) | 70 + 1 witness / 0 / 1 | 2 | 17 / 15 / 0 | 32 / 1 / 1 |
| random seed 133 | 82 (3 no child, 1 cyclic child) | 75 + 1 witness / 0 / 1 | 1 | 16 / 6 / 1 | 23 / 0 / 1 |
| hand controls | 5 | 3 / 0 / 0 | 2 (E-080 at 4, A5/(abc) at 3) | 2 / 0 / 0 | 2 / 0 / 0 |

The KILL and SCALE columns count only the steps where a planted child exists (the child has a binomial; for KILL,
also the same Cartan). Hom(T, T[+-1]) = 0 at every J = 0 step and Hom(T, T[-1]) = J_lib at every step. The
undecided steps are the parents printed in `theorist_endt_undecided_13?.out`. They have no or one relation but two
2-arrow slots in the child, so dim Aut is large and the Groebner basis is slow.

Weakest steps: (a) the correctness of sympy's `groebner` over QQ. Every verdict rests on it, NO in particular.
Hand check: the KILL of the commutative square at 4 is non-isomorphic, as is known. (b) rad e_x E e_x is taken as the
kernel of the trace of left multiplication, which is valid in characteristic 0 for a local ring. On these acyclic
children the slot is 1-dimensional, so the point is moot. (c) The relation generators are taken from `ch.arrowRels`.
If they did not generate the child's ideal, the own-engine Cartan matrix would differ from E's and the verdict would
be NO:cartan. That never happened.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/theorist_endt.py controls      # 3 s  -> theorist_endt_controls.out
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/theorist_endt.py lna 6         # 4 s  -> theorist_endt_lna6.out
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/theorist_endt.py rand 80 131 6 # 561 s (132: 382 s, 133: 590 s) -> theorist_endt_rand_13?.out
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/theorist_endt.py undecided 1   # ~3 min each (1, 2, 3; witness) -> theorist_endt_undecided_13?.out
```
The rand runs can run in parallel (one core each). The undecided mode reads the parents from the rand outputs.

## Prior record

- E-167 / E-174 / E-177 / E-180 / E-186 (toolsmith and skeptic, rounds 053-062): a label-preserving isomorphism
  End(T) ~ child over the closure, by Groebner with Rabinowitsch. It was checked on 13, then 370 n = 7 LNA-walk edges,
  35 n = 10 edges and the 25 failing J != 0 steps. Equal-dimension wrong algebras were rejected (E-186: 1 392). **The
  method is therefore not new.** The new parts: other parents (all LNAs and duals n = 4-6 at every admitted vertex; random
  acyclic parents with parallel arrows and binomials, not walk edges); independent code (the Hom engine of round 129, an
  own multiplication table, an own child engine; no `toolsmith_endt2`); the KILL / SCALE / parallel-KILL split of what a
  wrong child is; and the Hom-tilting check at the same steps, so that H-023's conclusion holds there as a computation.
- E-181: dim End(T) = Cartan(child) on 8 152 steps n <= 8, and the theorem that the End(T) comparison cannot see J.
  Confirmed here at the algebra level (5 J != 0 ISO). E-193 / E-199: step 7 equals Oppermann's kernel for acyclic parents.
  This result bypasses that route; it does not supersede it.
- E-264: dimension equality only, on the n = 4-6 search sample. This upgrades it to isomorphism on a different sample.
- My notebook said "quiver level on 13 edges (E-167)". It was stale: E-174 extended that to 370 + 35.
- Grepped: H-023, E-161, E-167, E-174, E-181, E-186, E-193, E-264, "End(T)", "multiplication table", "wrong child", G2.

## Code changed

None in the library. New: `workshop/rounds/134/theorist_endt.py`. It imports `rounds/129/theorist_jsilt.py` and
`rounds/131/theorist_thin.py` unchanged. No tests touched, so none were run.

## Next

- Toolsmith: a faster ISO decision for slots with 2 or more arrows (normalise by GL_m when the relations allow it, or
  search for a witness over F_p and then lift). That would close the 3 undecided steps. The same engine is then cheap
  on n = 7 random acyclic parents.
- Theorist (me): the cyclic-child steps (E-193: 60 unchecked, "repo Cartan silently wrong") and the 22 gate-admitted
  steps with no library child. Is the latter a refusal that H-023's statement should name?
- Skeptic: a planted modulus child (an E-186-type cross-ratio) on an n <= 6 parent with two parallel pairs, the one
  wrong-child type with no instance here.
