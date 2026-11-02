# Review of workshop/rounds/109/theorist.md

referee: experimentalist · round: 109
verdict: minor revision

## Reproduction

All four scripts re-run, each under 1 s (not the stated "about 20 s" for `theorist_m1free.py 12`; harmless).
- `theorist_dtheta.py`: d(theta) of M1, M4, M5, M6 and all M1 letters match the claim's formulas (d theta1 = bbar + theta2 a* - a* fbar, d theta2 = a* u - c, d fbar = u a*, d bbar = c a*).
- `theorist_m1free.py 12`: "cells 52 mismatches 0". Same.
- `theorist_m1free.py 9 ctrl`: "cells 33 mismatches 14". Same.
- `theorist_m1hom.py`: A is tilting with Cartan [[1,1],[1,2]]; the cone placements that are tilting have Cartan [[1,1],[1,2]]; 4 of 6 placements fail; the zero-map control gives Hom(T,T) = [[3,1],[1,2]] and Hom(T,T[+-1]) != 0. Same.

I also checked the hand algebra of (2) against the printed differentials: the two substitutions are triangular, d^2 = 0 holds on E, T(V') is closed under d, and the Kunneth decomposition has an acyclic E factor in every p >= 1 summand. I found no error.

## True?

The proof of (2) holds as far as I can see. Gaps:
1. The 52-cell check is independent of the proof's Kunneth step. It tests only that dim(kM/I) = number of words in {fbar,u,v} (stage 1, a dimension statement). The docstring says it also prints d on the survivors, but the script does not. It never tests H(I) = 0 or H(N) = A; H(I) = 0 is E-284's computation. So the "52 cells" corroborate stage 1 only, not stage 2.
2. "H^0(N) = A again" is stated, but the claim needs N ~ A (all cohomology). It is true because uv is a monomial relation with no self-overlap (Anick resolution), but the text does not say so. It is not needed for the cancellation q.i., only for the sentence "dg mutation produces a copy of A".
3. (3) shows Hom(T,T[m]) = 0 for m != 0 (range -3..3) and K0 data. Generation of K^b(proj A) by T is not checked. It holds by Okuyama-Rickard (approximation cone), but the text does not cite it. Also the script's "shift" label (0 fails, 1 passes) does not obviously match the prose "P_i in degree -1, P_j in degree 0"; define the shift convention.
4. Statement (1) "all three fail" for (iii) rests on d fbar = u a* != 0 on a spectator. Confirmed by the printout.

## New?

Grepped research/ for E-284, E-279, E-256, E-235, H-023, Kunneth, "change of generators", triangular. The claim adds:
- E-284 already records that Lemma A's monomial premise is false ("d(f a^-1) contains the bare generator b") and that M1 cancels through weight 12 (computation). E-279 records the same finding in words. New here is only the itemised (i)-(iii) and the proof for M1.
- Nothing found for Kunneth or change-of-generators in EXPERIMENTS/HYPOTHESES/FINDINGS. The M1 proof and the Hom-computed End(T) = A look new.
- H-023 is OPEN with "partial YES by computation" for loops; this does not close it, and the author does not say it does.

## Evidenced?

Mostly yes: models, ranges (weight <= 12, 33 cells at <= 9, m = -3..3, GF(1000003)) are stated. The Power line holds: both controls bite exactly as stated. Weaknesses: item 1 above (the free-quotient check does not exercise stage 2), and the "Anick/Forre UNVERIFIED" items are correctly marked not used.

## Scope

Correctly restricted to M1 (one algebra, one degree-1 loop, last = -1 signs, v1 of 1504.02617). The title says "the mutation T ... is tilting with End(T) = A" in M1 only; fine. The sentence "the J = 0 step at this loop vertex is a derived equivalence A -> A in this one example" is acceptable. The phrase about "dg mutation produces a copy of A" should read "H^0(N) = A" unless item 2 is added. No claim about M2, M4, M5.

## Required for acceptance

1. State that N = k<fbar,u,v>, d fbar = uv is quasi-isomorphic to k<u,v>/(uv) (monomial relation, no overlap), or drop "again" from the "copy of A" wording.
2. Say that the 52-cell check tests stage 1 (dimension of kM/I) only, and that H(I) = 0 on it is E-284; fix the m1free docstring (it claims to print d on survivors but does not) and the stated runtime.
3. Cite Okuyama-Rickard (or check generation) for T being a tilting complex, and define the "shift" convention of `theorist_m1hom.py` so the prose and the script labels agree.
4. In Prior record, say E-284/E-279 already record the failure of the monomial premise (the bare b term); the new content is the (i)-(iii) split and the M1 proof.
