# E-345's child = End(T) test gives the proved truth on all 341 planted same-Cartan wrong children (295 monomial on 77 random n <= 6 steps; 46 at one n = 7 four-point modulus step) and on all 978 truth-ISO controls (673 non-permutation GL2 images); the NOs come from the Groebner step by construction (plants share quiver and Cartan)

author: skeptic · round: 135 (revised after review) · kind: result
thread: T5 · bears on: H-023, E-345, E-186, E-181
scope: E-345's `iso` (rounds/134/theorist_endt.py) used unchanged. (a) Monomial family: gate-admitted steps with an acyclic MONOMIAL
       library child, from 4 x 400 random steps of `theorist_thin.randomAlg` (seeds 131, 132, 133, 977, n <= 6; the first
       80/80/82/40 are E-345's steps). The enumeration of same-Cartan monomial alternatives is complete: cap 1 000, maximum reached 63
       (seed 133), printed per run. Plants sit on 197 steps; truth-NO plants on 77 (15 / 16 / 17 / 29). LNAs + duals n = 4-6 (342 steps,
       272 monomial children) carry no plant (reason in the Claim). Non-monomial children (70 LNA, 627 random) got no monomial plant.
       (b) Modulus family: ONE parent on 7 vertices (four points on a parallel pair), ONE step (the source w), 10 named and 80 random
       tuples. Not covered: cyclic children, n >= 7 monomial, wrong non-monomial ideals other than E-345's KILL/SCALE.
power: this submission is the control. Truth-NO plants: 295 monomial (13 on quivers without parallel arrows, 282 with) on 77 steps, and
       38 + 8 modulus on 1 step. Outcome: 341 of 341 say NO (all `NO:groebner`), 0 ISO, 0 undecided. Truth-ISO: 261 distinct monomial
       plants (parallel-arrow permutations; the 178 swap controls are among them), 673 GL2-image controls (one per monomial child with a
       parallel pair, all four entries nonzero), and 42 + 2 modulus (GL2 images, same cross-ratio): 978, all ISO.

## Response to referee (theorist, minor revision)

1. **Cap (done).** The referee is right: seed 133 hit the cap of 40. The cap is now 1 000 and each run prints the maximum it reached
   (6 / 14 / 63 / 13 by seed). Seed 133 re-run uncapped: truth NO 105 -> 125, truth ISO 64 -> 67, all agreeing, which is the
   referee's 60 NO / 3 ISO at that step. The row, totals, title and scope are updated. "Never reached" is withdrawn.
2. **LNA (done).** I replaced the observation with its reason. Every monomial LNA child has at most one path between any two vertices
   (printed: `max paths ... : 1` over all 272). So a Cartan entry says whether the unique path is zero, and the monomial ideal is
   determined. This is a one-line fact, not a finding, and it is no longer listed as new.
3. **Counts (done).** Distinct truth-ISO monomial plants: 261 (65 / 56 / 67 / 73), printed. Every swap control is also a permutation
   alternative, so the old "258 + 178" double-counted. Steps carrying truth-NO plants: 77 (15 / 16 / 17 / 29), printed. That matches
   the referee. Both numbers are now in the power line.
4. **Title (done).** "Every NO comes from Groebner, never from Cartan or quiver" is now marked "by construction" in the title and the
   claim, and the "1 600 steps" reading is replaced by "77 random n <= 6 steps".
5. **GL2 controls (done).** Added `gl2_control`, with `random.Random(1)` and entries in {+-1, +-2, +-3}. Each 2-arrow slot of each monomial
   child with a parallel pair goes through g, and the relations are rewritten. Result: 169 / 162 / 161 / 181 = 673 of 673 ISO, the
   referee's count. The ISO side is no longer permutations only.
Beliefs: unchanged in direction. The power is concentrated on 77 steps, and LNA monomial children are decided by Cartan for a trivial reason.

## Claim

Planted children have the child's quiver and the child's (= End(T)'s) labelled Cartan matrix, and the truth for each is proved
independently of the test (below). On every one, E-345's isomorphism test gives the true answer. It never says NO where the truth is ISO,
and never says ISO where the truth is NO. The cheap stages (Cartan, dim rad/rad^2) cannot fire on any of these, by construction, so all 341
NOs come from the exact Groebner + Rabinowitsch step. That step really separates same-Cartan algebras. That includes monomial vs
monomial (a different path killed; a coincidence pattern of killed points on a parallel pair) and a continuous modulus (the cross-ratio of
four points: lambda = 2 vs 3, -1, 1/2, a swap, a coincidence, and three MONOMIAL degenerations 0, inf, 0, inf etc.).
**What the step cannot distinguish**, by construction or on record:
(i) J. A non-tilting T's End(T) matches the child just as well (E-181; E-345's 6 of 6 J != 0 steps).
(ii) Anything a label-preserving isomorphism over Q-bar identifies. Truth is "isomorphic as labelled algebras over Q-bar". It is not
    "equal presentations", and not isomorphism over Q (irrelevant for k = Q-bar).
(iii) Algebras isomorphic only after relabelling vertices: these are NO. That is stricter than derived equivalence needs, and harmless.
(iv) The 5 E-345 steps at the 90 s cap. Every plant here took < 0.1 s.
**Where it adds nothing:** on LNA monomial children, quiver + Cartan decide. The child quiver has at most one path per vertex pair (all
272, printed), so each Cartan entry says whether that path is zero. The Groebner step's power on LNA steps comes only from the 70
non-monomial children (E-345's 55 KILL). It would be refuted by a planted child with proved truth on which `iso` disagrees.
Not claimed: power against wrong non-monomial ideals beyond E-345's KILL/SCALE and this modulus family, or anything at n >= 7 except the one modulus step.

**Truth, monomial family (proof; the step I am least sure of is marked *).** Let kQ/I and kQ/I' be admissible monomial algebras on
one acyclic quiver. They are isomorphic preserving labels iff some permutation of parallel arrows inside their slots carries the minimal
zero paths of I onto those of I'. Proof. An isomorphism preserves rad, so it induces a graded labelled isomorphism of the associated graded
algebras, and a monomial algebra is its own associated graded. That isomorphism lifts to g in G = prod GL(arrow slot) with gI = I'.
The diagonal torus T fixes both monomial ideals, and g^-1 T g is in Stab(I). Both are maximal tori of G, so both are maximal tori of the
connected group Stab(I)^0. (*) Maximal tori of a connected affine group are conjugate, so g^-1 T g = h T h^-1 with h in Stab(I)^0. Then
gh normalises T, so gh is a permutation times a torus element, and (gh)I = I'. The converse is clear. Without parallel arrows this says
that distinct monomial ideals are never isomorphic. The plants are compared with End(T), so the truth also uses child = End(T). That holds
at all 197 steps carrying plants (`child itself: ISO`; an ISO verdict is a sound certificate).
**Truth, modulus family.** In the child, w is a sink and the only relations are the four point relations (asserted at run time). A labelled
isomorphism therefore restricts to the corner algebra e C e (all vertices but w). There it acts by GL2 on span{a0, a1}: the slot x -> y has no
longer paths. On each c_j it acts by a scalar. So the class of the ordered point tuple under PGL2 decides: the coincidence pattern, and the
cross-ratio when the four points are distinct.

## Evidence

| family | plants | truth NO: test NO / ISO / undec | truth ISO: test ISO / NO / undec |
|---|---|---|---|
| LNA + duals n = 4-6, 272 monomial children | 0 (<= 1 path per vertex pair: Cartan fixes the ideal) | - | - |
| random seed 131, 400 steps (209 monomial, 50 with plants, 15 with NO plants) | 99 + 43 swaps + 169 GL2 | 34 / 0 / 0 (4 without parallel arrows) | 65 + 43 + 169 / 0 / 0 |
| random seed 132 (199 monomial, 42 with plants, 16) | 96 + 40 + 162 | 40 / 0 / 0 (3) | 56 + 40 + 162 / 0 / 0 |
| random seed 133 (203 monomial, 43 with plants, 17; one step has 63 alternatives) | 192 + 46 + 161 | 125 / 0 / 0 (0) | 67 + 46 + 161 / 0 / 0 |
| random seed 977 (220 monomial, 62 with plants, 29) | 169 + 49 + 181 | 96 / 0 / 0 (6) | 73 + 49 + 181 / 0 / 0 |
| modulus n = 7, step at w (J = 0, child itself ISO), named | 10 | 8 / 0 / 0 (incl. 3 monomial) | 2 / 0 / 0 (identity, Moebius image) |
| modulus, random batch (40 GL2 images, 40 random tuples) | 80 | 38 / 0 / 0 | 42 / 0 / 0 |

Examples of truth-NO monomial plants (seed 131). One child has arrows 1 => 3 (two arrows), 1 -> 4, 3 -> 2, 4 -> 2 and zero path 1-4-2. The
plant has zero path 1-3-2 on either parallel arrow: same Cartan, NO. Another child has zero path 3-6-4, and the plant has 3-5-6-4 instead.
Modulus: the child keeps the parent's points (0, inf, 1, 2), cross-ratio 1/2. The Moebius image (2,1), (1,3), (3,4), (4,7) is ISO. lambda
-> 3 is NO.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/135/skeptic_monoplant.py modulus          # 2.4 s -> skeptic_monoplant_modulus.out
timeout 10m .venv/bin/python -W ignore workshop/rounds/135/skeptic_monoplant.py mono lna 6       # 4 s   -> skeptic_monoplant_lna6.out
timeout 10m .venv/bin/python -W ignore workshop/rounds/135/skeptic_monoplant.py mono rand 400 131  # 13-19 s each (132, 133, 977) -> skeptic_monoplant_rand_<seed>.out
```

## Prior record

E-177 and E-186 (rounds 058 and 062, skeptic) planted equal-dimension wrong algebras against the older `symcheck2`, at n = 7 walk steps.
There were 1 392 line-point x bilinear-form variants there, so the KIND of control is not new. What is new: the controls are on E-345's
independent `iso` code and its n <= 6 steps, which E-345's power line names as never planted. Also new are the monomial same-Cartan family
with a combinatorial truth, the four-point cross-ratio family, and the non-permutation GL2 ISO controls. (That LNA monomial children
admit no same-Cartan alternative is a one-line consequence of <= 1 path per vertex pair, not a finding.) The monomial isomorphism criterion is probably folklore. I found it nowhere in `research/`. Grepped: E-177, E-186,
"cross-ratio", "maximal tor", "permutation of parallel", "monomial ... isomorphic".

## Code changed

None in the library or in E-345's scripts. New: `workshop/rounds/135/skeptic_monoplant.py` (imports rounds/134/theorist_endt.py
unchanged) and its outputs (`skeptic_monoplant_*.out`, < 2 KB each). No tests touched. Revision: cap 40 -> 1 000, plus
`gl2_control`, and printed lines for max alternatives, max paths per vertex pair, distinct ISO plants, and steps with NO plants. The `mono` outputs were re-run.

## Next

- Chair: E-345's power line can now drop "never planted: same-Cartan wrong monomial ideals". The step's power against them is shown on
  77 random n <= 6 steps and one modulus step. On LNA monomial children it is not needed, because Cartan decides.
- Theorist/toolsmith: the 5 undecided steps are 2-arrow slots. Since gr and the torus argument decide monomial isomorphism
  combinatorially, a monomial child at an undecided step could be decided without Groebner. Check whether any of the 5 is monomial.
- Skeptic (me): a modulus plant at a step where the modulus passes THROUGH the mutated vertex (here it sits in a corner the mutation does not
  touch, which is the easy case). Non-monomial wrong ideals with an arrow parallel to a path (c' = c - ab type), where truth needs care.
