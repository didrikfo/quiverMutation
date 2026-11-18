# Review of workshop/rounds/131/theorist.md

referee: experimentalist · round: 131
verdict: minor revision

## Reproduction

All six commands re-run with `timeout 10m`.

- `theorist_thin.py planted`, `sweep 600 131`, `cyclic 600 131`: 3.5 s together. Every line matches the submission: E-080 square (gate True, dmax 2, J 1, lib (True, False)); ab = 0 control; monomial square; E-318 parent v = 3 d = 2 not tilting, v = 1, 2, 4, 6, 8, 9 tilting. Both sweep tables are the same count for count. Library mismatches 0, L1 violations 0. Cyclic (1 739, 1 739, 380).
- `theorist_thincyc.py` on the 13 pairs: 71 s in total. Per-pair calls are 321, 943, 361, 181, 212, 181, 421, 698, 617, 435, 943 at n = 10 and 181, 576 at n = 11, which sums to 6 070. Every pair has d2 0, dmax 1, refused 0, cyclic 0. This matches.
- One case further, with a new seed and about 3 to 5 times the size (6.5 s). `sweep 2000 7`: the mixed sample has (2,T,T) 92 and no (<=1,T,T) row. The monomial sample has (2,T,T) 0. No row has gate F with J = 0. Library mismatches 0, L1 violations 0. `cyclic 3000 7`: (8 823, 8 823, 1 894). Path-simplicity held at all 8 823 thin vertices.
- Record cross-check for (c): I tested E-329's counterexample to the library gate (`toolsmith_cycgate.py` bag: a:4>5, b:5>4, g:5>1, x1:1>2, x2:2>3, relations ab and bagx1, at i = 1). By hand, e_5 A e_1 = span{g, bag}, so d_5 = 2. The vertex is thick, and the refusing path bag is the non-simple path that (c) predicts there. This is consistent with (c), and it is the w != v case of the proof on a recorded instance.

## True?

(a) is correct line by line, and so is (b) (trivial direction).

(c): both cases are correct. The scalar argument needs only C^m in I. The claim states R nilpotent mod I as a hypothesis, so the result stands. The "weakest step" remark is false as written, though. Finite-dimensional A together with I in R^2 does not give R^m in I. Counterexample: one vertex, one loop x, I = (x^2 - x^3), which lies in R^2. Then A = k[x]/(x^2(1 - x)) has dim 3, and x^m != 0 for every m because the factor k[x]/(1 - x) sees x = 1. A loopless analogue on a 2-cycle (I = (ab - abab)) should behave the same way. Nilpotence must stay a separate hypothesis (I admissible). It does not follow from finite dimension.

(d) is correct. Fix a lift of x supported on paths that are nonzero in A. Monomiality of e_i I e_t(beta) then puts each p beta in I, and the gate refuses p. Unlike (a), this per-row statement is not computed: neither script tests whether some e_i I e_t(beta) is non-monomial at the 29 (or 92) (2,T,T) rows. Only the monomial-sample corollary (0 rows) is computed.

## New?

- (a) is E-128 (L1), printed there with "J_i != 0 forces d_i >= 2". (b) is the definition of the gate read backwards, and it sits inside E-128's "rejects v iff" sentence. J = 0 being equivalent to tilting is E-130 / E-171.
- The monomial half of (d) is E-208 / E-310. The d = 2 binomial shape is E-080 / E-318 / E-131.
- (c) is new. E-130 states the caveat "the repo gate enumerates simple paths only", E-128 / E-171 is cited for the `allPathsBetween` caveat (E-213-era line), and E-329 records library-gate errors with a cycle upstream. I found no entry that ties exactness of the simple-path gate to thinness (grepped `vertex-simple`, `simple paths only`, `non-simple path`, `nilpoten`, `thin`).
- The acyclicity count on 13 F-041 pairs is new data, but it is small. It explains E-337. It does not extend it.

## Evidenced?

Yes. The proofs carry the claim, and the computations are exact (own row-reduction over Q). Power is shown: a planted failure fires for each dropped hypothesis, and I re-ran them. The cyclic sample is monomial with rad^4 = 0. That makes it a combinatorial check of (c), with no library comparison on cyclic quivers, which the submission disclaims. The 380 / 1 894 thick non-simple counts show the test is not vacuous.

## Scope

The title matches what is proved. The F-041 acyclicity statement covers 13 of 1 817 pairs, and the submission says so. "Any field": the proofs use no field property, but the computations are over Q only. That is acceptable for a proof.

## Required for acceptance

1. Correct the "weakest step" sentence in the proof of (c). Finite dimension plus I in R^2 does not give R^m in I (counterexample k[x]/(x^2 - x^3)). State that (c) needs R nilpotent mod I (I admissible) as its own hypothesis.
2. Either compute (d) per row or narrow its wording. For each (2,T,T) row of `sweep 600 131`, print whether some out-arrow beta has e_i I e_t(beta) not spanned by paths, and report the count (should be 29 of 29). Otherwise say that only the monomial corollary is computed.
3. (Optional, one line.) Add E-329's bag instance as a record check of (c): it is thick (d_5 = 2), and its refusing path bag is non-simple.
