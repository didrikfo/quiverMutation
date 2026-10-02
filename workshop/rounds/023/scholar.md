# "Reject iff long-sided square" splits: long square => reject needs only minimality; reject => long square needs out-degree 1 (a hypothesis that is not a theorem) and a choice of generators

author: scholar · round: 023 · kind: proof · thread: T5 · bears on: H-015, E-066, E-097, E-100

## Claim

Notation as E-097: v the mutated vertex, alpha the arrows out of v, J_i = {c in e_iAe_v : c alpha = 0 in A for every alpha}, I the ideal.
Facts used (E-097): tiltingPlus False <=> some g_i not injective <=> J != 0 (as a quotient of J by I, i.e. c not in I). The gate (`procedure.isMutable`,
read from its code) refuses exactly when some *single path* p not in I has p alpha in I for all alpha. So "reject" (gate admits, tiltingPlus False) means:
J has a nonzero class, and no nonzero class is a single path. Step 7 enters only to say that the code's rewrite fails the Cartan test exactly then
(E-097 (b), step-7 completeness); the shape statement itself below does not use step 7.

(1) Long square => reject. Hypothesis: v has exactly one outgoing arrow alpha, and I has a *minimal* generator r = c alpha with c a combination of >= 2
paths into v (c = sum of the paths of the relation with the last arrow alpha removed). Then c alpha in I, and c not in I, because otherwise r = c alpha
lies in I alpha and is not a minimal generator. So J != 0, tiltingPlus is False; the gate admits because c is not a single path (and no path is killed,
else the gate refuses and this is a gate refusal, not a step-7 rejection). Needs: minimality of r; out-degree 1 is NOT needed for this direction
(with several arrows one needs c killed by each, see D).
(2) Reject => long square. Given c in J, c not in I, c alpha in I. With out-degree 1 this says: some element of I of the form c alpha has c not in I, i.e. I is not
(I : alpha) = I. Hypotheses needed: (a) out-degree of v is 1; (b) "a relation in the presentation" is read up to a change of generators (G below: the
presentation {P1 - X, P2 - X} has c = P1 - P2 but no generator whose paths all end in alpha). (a) is an empirical regularity of the algebras reached, not
a consequence of the gate or step 7: case D (below) is rejected with out-degree 2.
Not claimed: (a) holds on all gate-admitted algebras; the hypothesis "step 7 complete" (E-097 (b)) is the only unproved link between this and the code's behaviour.

Why E-078's length-3 square is fine (my notebook's question): then the relation lies at v itself, c in I (c alpha = 0 trivially, class 0); a relation reaching
past e has c alpha not in I. The kernel needs the relation to end exactly at the arrow(s) out of v. That is what "long-sided" means.
Why no path kernel survives on walks: the gate removes the monomial case (F). So the commutativity element is forced: a >= 2-term relation, as in E-066.

## Evidence

Hand-built cases (a=1,b=2,c=3 square into v=4; e=5,f=6; X path 1->7->5), `scholar_longsquare.py --hand`:
| case | gate | out(v) | dim J | tiltingPlus | long square in rels |
|---|---|---|---|---|---|
| A relation abve = acve | admit | 1 | 1 | False | yes |
| B relation abv = acv (short) | admit | 1 | 0 | True | no |
| C relation abvef = acvef (past e) | admit | 1 | 0 | True | no |
| D commutes into e and into f | admit | 2 | 1 | False | no (out 2) |
| E commutes into e only | admit | 2 | 0 | True | no |
| F zero relation abve = 0 | refuse | 1 | 1 (single path) | False | no |
| G {abve = a x e, acve = a x e} | admit | 1 | 1 | False | no (not all paths end in v,e) |
D is a reject with out-degree 2, so (2)(a) is a hypothesis; G is a reject that the literal `alg.rels` test misses (a presentation artefact); F is where the gate, not step 7, speaks.

Guided walks (class 0, gate-admitted steps, 300 s cap each; counts are steps, not distinct parents; table `scholar_longsquare_n{5,6,7}.txt`):
| n | steps | J != 0 | of those: out-degree 1 / long square on alg.rels | J = 0 with long square | J != 0 with out-degree 2 |
|---|---|---|---|---|---|
| 5 (closed, 6 240 algebras) | 16 620 | 0 | - | 0 | 0 |
| 6 (46 858) | 101 045 | 1 139 | 1 139 / 1 139 | 0 | 0 |
| 7 (35 094) | 66 049 | 162 | 162 / 162 | 0 | 0 |
So on the walks "J != 0 <=> out-degree 1 and long square on alg.rels" holds step by step (consistent with E-100, which counted the same on 479 761 steps); the
walks reach no D or G type algebra. Out-degree 2 occurs in 33 309 + 17 933 admitted steps, all with J = 0.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/023/scholar_longsquare.py --hand                       # 2 s
timeout 10m .venv/bin/python workshop/rounds/023/scholar_longsquare.py 5 --class 0 --budget-sec 300   # 26 s
timeout 10m .venv/bin/python workshop/rounds/023/scholar_longsquare.py 6 --class 0 --budget-sec 300   # 300 s (cap-dependent)
timeout 10m .venv/bin/python workshop/rounds/023/scholar_longsquare.py 7 --class 0 --budget-sec 300   # 300 s (cap-dependent)
```

## Prior record

E-066 located the commutativity element; E-097 gave kernel/cokernel and the long-sided shape; E-100 gave the walk control (no tilting step has it, all 2 104 rejecting
do) and left "is it a theorem" open. New here: the two directions with their hypotheses, the gate's single-path/combination split as the reason the kernel
is a >= 2-term element, the "short level is fine" one-liner, and the out-degree-2 example D showing the shape is not forced. Nothing here is in RETRACTIONS.
The step that is still not proved is E-097's step-7 completeness; the (2) direction is a restatement of J != 0 once out-degree 1 and generators are fixed.

## Code changed

None in `quivermutation/`. New `workshop/rounds/023/scholar_longsquare.py`; no tests touched.

## Next

- toolsmith/theorist: is out-degree >= 2 with J != 0 reachable from an LNA by guarded steps (D-like)? A bounded search at n = 8, 9 for gate-admitted J != 0 with out-degree >= 2 would decide whether "long square" is a theorem on the walk class or only a pattern.
- theorist: referee (1) and the step-7 completeness link; test case D through the actual rewrite (Cartan congruence should fail).
- chair: wording for E-084/E-095/E-100: "long-sided square" is a shape of the walk class, with the iff holding intrinsically as J != 0.
