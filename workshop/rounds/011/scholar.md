# A 5-vertex commutative square with a length-4 relation is a gate-admitted, non-tilting mutation: the second such case, and it is tiny

author: scholar · round: 011 · kind: result
thread: T5 · bears on: H-015, E-032, E-055, E-057, E-066, STEERING round 002 q1 (`isTilting`)

## Claim

On the hand-built algebra A5 (vertices a,b,c,d,e; arrows a>b, a>c, b>d, c>d, d>e; one relation
abde = acde, so e_a A e_d has dimension 2 and the composite with d>e kills the commutativity
element), the gate `mutationIsPossibleAtVertex` admits d, `tiltingPlus` (the one map of AI 2.32(b) =
Ladkani 2.3(c)) returns False, and the repo's rewrite fails the Cartan congruence
Cartan(child) = R C R^T. So the gate-admitted rejection of E-032 step 7 is not special to n = 10:
the smallest instance has n = 5. The same happens at n = 6, 7 with the instance padded by a chain
before a or after e (6 of 6 padded placements at n = 5..7; 1 gate-admitted rejection per algebra,
always at d). It does NOT claim the A5 algebra is reachable from an LNA by guarded mutation (not
tested), nor that the E-032 parent is the same case up to padding (it is the same shape, not
checked as an isomorphism).

CHZ 2509.12983: **still UNVERIFIED.** `curl https://arxiv.org/pdf/2509.12983` and `/abs/` both
return "CONNECT tunnel failed, response 403" (organization policy; the proxy status names
arxiv.org:443 connect_rejected). I did not try other hosts or disable anything.

## Evidence

Script builds, for n = 5..7, a chain prefix (pre) + square + d>e + chain suffix (post), with three
relation kinds, and at every vertex with an outgoing arrow records (gate, tiltingPlus); where the gate
admits, also the Cartan congruence of the actual rewrite.

| kind | relation | n = 5, 6, 7 instances | at d: gate / tiltingPlus | other vertices |
|---|---|---|---|---|
| long | abde = acde | 1 + 2 + 3 = 6 | True / **False** (cong False) | all (True, True), cong True |
| short | abd = acd (true square) | 6 | True / True | all (True, True) |
| zero | abde = 0, acde = 0 | 6 | False / False | all (True, True) |

Totals over the 18 algebras: gate-admitted with tiltingPlus False: 6, all `long`, all at d, all
also fail the congruence (so `tiltingPlus` and the congruence agree on every admitted step here, 0
disagreements). Gate-refused with tiltingPlus True: 0. The `zero` row is the monomial control: the
gate itself refuses d there, as E-057's control says. The `short` row shows it is the relation
length, not the square, that matters: when the commutativity is at length 3 the map on e_a A e_d is
injective (dim 1) and the step is tilting.

Mechanism, one line: p |-> (p d>e) kills the 2-dimensional e_a A e_d to dimension 1 exactly when
the relation passes through the outgoing arrow. This is E-066's kernel element, now at n = 5.

CHZ wording, argued but not read: Prop 3.5 uses tail-maximal paths (soc P_i). In
`research/literature/2509.12983` Cor 3.6 is stated path-wise ("nonzero path p prolongs to a nonzero
path"). For a non-monomial I "p not in I" is meaningless for a path in a quotient where
abde = acde (each of the two paths is nonzero, their difference is zero), and the E-066 parent
passes the path-wise wording while failing Prop 3.5. That is a reason to expect a monomial (or at
least "paths are a basis") hypothesis; it does not show the paper's text lacks one. The note's
header says it was read from the arXiv PDF on 2026-09-19, so someone with the PDF can settle it in
one look at the statement of Cor 3.6.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/011/scholar_square.py     # under 10 s
```

## Prior record

E-066 listed "n = 10 is the first size with a commutative square into a vertex with one outgoing
arrow" as untested and asked for a smaller instance: this answers it (n = 5). E-055/E-057 found no
gate-admitted rejection at n <= 7 because their parents were LNAs and their relation duals walked
from LNAs, and no such algebra was reached; that is consistent with this result, since A5 is not
shown reachable. Whether any gate-admitted rejection occurs on a path from an LNA at n <= 9 is
exactly the open audit (Menu 4). The one-map identity itself (AI = Ladkani = `tiltingPlus`) is still
the author's statement, now with 18 more algebras where `tiltingPlus` and the Cartan congruence of
the true rewrite agree, which is a check of the code against the actual rewrite, not of AI
against Ladkani.

## Code changed

None (one new script, `workshop/rounds/011/scholar_square.py`, reusing `scholar_h015.py` helpers).
No tests run, since no library file was touched.

## Next

- Chair: STEERING round 002 q1 asks for "a gate-admitted rejection" before promoting `isTilting`;
  A5 is one, though hand-built. Decide whether a hand-built non-LNA counts or whether it must come
  from a guarded walk. Either way a unit test with A5 (gate True, tiltingPlus False) is cheap.
- Toolsmith: is A5 (or its padded forms) reachable from an LNA by gate-admitted steps at n <= 9?
  If yes, the walks of E-055 at larger n should hit it; if no, say why the gate-admitted rejection
  needs the tree-with-square shape.
- Anyone with PDF access: read the statement of CHZ Cor 3.6 for "monomial" and update the flag in
  `research/literature/2509.12983`.
- Theorist: derive the one-map identity (AI 2.32(b) vs Ladkani 2.3(c)); independent End(T) Cartan
  for A5 at d is a small hand computation.
