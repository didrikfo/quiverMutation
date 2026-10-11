# The literature confirms, but does not independently test, E-032 step 7: AI Thm 2.32(b) is the same linear map as Ladkani 2.3(c), and it fails by a hand-checkable commutativity element

author: scholar · round: 006 · kind: negative · (no new data; reading plus two small scripts)
thread: T5 · bears on: H-015, H-010, E-032, E-057, E-059

## Claim

The rejection at E-032 step 7 is correct: the parent there has a nonzero element
c = [8,6,4] + [8,10,4] of e_8 A e_4 (paths from 8 to 4) with c * (4->9) = 0, which is exactly
the failure of injectivity in Aihara-Iyama Thm 2.32(b) for D = sum of P_j, j != 4, g built from the arrows at 4.
So the right mutation at 4 is silting, not tilting, and its End is not guaranteed derived
equivalent (the key does move, E-032/F-038). **It does not claim** that the literature
gives a check independent of the code: 2.32(b), Ladkani 2.3(c) and `tiltingPlus` are one
map (p |-> (p beta)_beta over arrows beta out of the vertex), so agreement is not corroboration
of the implementation. The literature does not decide H-015 (guard sufficiency): both
criteria are per-step exact, which is what E-057 already uses.

## Evidence

1. Replayed `[4,6,4,6,9,4,4,6]` from the relation dual of `03033030` (`scholar_step7.py`). Parent at step 7
   (10 vertices): arrows 1>2>3>5>7>8, 8>6, 8>10, 6>4, 10>4, 4>9; relations include the
   commutativity `8,6,4,9 + 8,10,4,9 = 0`. The only arrow out of 4 is 4>9. Over all i != 4 the map on e_iAe_4
   has full rank except i = 8 (dim 2, rank 1). Kernel = c above. Step 7 is the first step on this
   path where a commutativity relation reaches through the vertex; steps 1-6 and 8 pass.
2. Why the gate admits it: the gate (2112.08129 quiver criterion) is monomial-blind. No relation
   *ends* at 4's head in a way it sees; the obstruction is a non-monomial relation, a socle
   element that is a sum of paths. That is the mechanism, and it is why E-059's non-monomial
   parents at n <= 7 never showed it: n = 10 is the first size with a commutative square
   feeding a vertex that still has one outgoing arrow.
3. Side matters (`scholar_sides.py`): the mirror test (arrows into k) is False at steps 1, 3, 6 where the key
   holds and the right test is True, and True at step 7. So the repo performs the *right* mutation
   mu^-_{P_k}(A) and tiltingPlus tests the right side; the test on the wrong side would give
   spurious results both ways. (Convention of "right" is inferred from this, not from reading code.)
4. CHZ (2509.12983) Cor 3.6 with |S| = 1: read *only via the repo summary* (arXiv blocked by the
   proxy this round). As summarised, the path-wise wording ("every nonzero path prolongs to a nonzero path")
   would PASS at step 7: p1 = [8,6,4] and p1*(4>9) are each nonzero in A; only the sum dies. The
   underlying Prop 3.5 (Phi^+ via supp soc P_i) FAILS: c is a socle element supported at 4, so
   4 in Phi^+({8}), violating closure of S^c. So the path-wise specialisation is valid for monomial I
   (all LNAs, where the 2052-algebra agreement was measured), and **not** for non-monomial parents.
   I have not seen the paper's hypothesis on I; UNVERIFIED. The summary file should not say
   "Cor 3.6 is an iff for kQ/I" without "monomial" until someone reads the PDF.
5. CHZ does not decide either way for the search's completeness: Cor 3.20 is for two-term complexes
   over the current algebra and, as the summary says, does not compose along a chain.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/006/scholar_step7.py   # seconds (needs no RD env; relation dual built inside)
timeout 10m .venv/bin/python workshop/rounds/006/scholar_sides.py   # seconds
```

## Prior record

E-057 and E-059 (Ladkani 2.3(c) agrees with gate; step 7 the only rejection), E-032 part 5 and
F-038 (step 7 moves the key), literature notes 1009.3370 (Thm 2.32, "implement this") and
2509.12983. New: the explicit witness, the AI = Ladkani = `tiltingPlus` identification, the
side check, and the monomial caveat on CHZ Cor 3.6. Not new: that the step is rejected.

## Code changed

None in `src/`. Two scripts in `workshop/rounds/006/` (import round-001 `scholar_h015.py`).

## Next

- Chair: promote a one-line caution to `literature/2509.12983` (Cor 3.6 path-wise form is monomial only; Prop 3.5 socle form is the general one), after someone reads the PDF (toolsmith/scholar with network).
- H-010 / T7 (theorist): step 7 fails at a commutativity element, so the "overlap reducible only at an end" proof
  from 2112.08129 must treat non-monomial relations; this is a concrete case to test it on.
- For genuinely independent corroboration: compute End of the two-term complex mu^-_{P_4}(A) directly
  (cone of D_4 -> P_4) and compare its Cartan matrix with the repo's child; toolsmith, small.
- Overnight audit n = 9/10 stays the only way to find a second gate-admitted rejection.
