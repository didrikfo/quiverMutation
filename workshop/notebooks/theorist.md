# Theorist notebook (rewritten round 109)

## What I believe now
- 109 (T5, loop case): Lemma A's hypotheses (i) single quadratic monomial d_1 theta, (ii) disjoint X1/Y/Theta, (iii) D = 0 on spectators all fail for a loop f at i. In M1 (A = kQ/(ab), d f = ab) d(f a^-1) = bbar + theta2 a* - a* fbar has a linear pivot bbar, so the Anick lemma is not needed: after the triangular substitution pi = d theta, T(V) = T(fbar,u,v) * T(E) with E contractible, and I = (E) is acyclic by Kunneth. Complete proof for M1 only (pivot-generator type); M2, M4, M5 have product or shared pivots and stay computational (E-284).
- Hom control: T = P_j + cone(P_i -> P_j) is tilting for that A (Hom(T,T[m]) = 0, m = -3..3), Cartan [[1,1],[1,2]], End(T) = A by hand = H^0(N). Planted errors (zero map, wrong shifts) say NO. So the loop step is a derived equivalence A -> A in this one example.
- 106 (still true): six hand-sized models cancel to weight 12; sign set 'last = -1'; no degree-0 loop is testable. 103: loop of degree >= 1 at i forces i onto an oriented cycle (E-279); Oppermann needs only no degree-0 loop at i. Anick numbers stay UNVERIFIED (secondary).
- Older: 101 cancellation is Prop. 7.6; 098 (T4) class A/B decided by first step (E-187/E-205); 093 conditional on omega-reduction.

## What I tried
- 109: printed d(theta) for M1, M4, M5, M6; hand elimination of pivots, checked by dim(kM/I) = words in the survivors (52 cells, 0 mismatches, planted control 14 mismatches); wrote a small Hom engine for complexes of projectives over kQ/(ab) (right modules, concatenation), 6 placements of the cone and a zero-map control.
- 106: loop filtration script; 103: loopdeg/cycsweep.

## Next
- Pivot-generator class: characterise it (loop f with d f = alpha beta, beta a single arrow) and prove for it in general; the shared-pivot case M4 by a change of basis in the theta span.
- Product-pivot cases (M2, M5): Theta-count filtration plus Lemma A, then discharge.
- Toolsmith: `allowCyclic=True` child of M1 vs End(T) = A; a Hom-computed End for a second loop algebra (e.g. relation loop of degree 2).
- T4: shape of the 14 width-5 middle algebras; 6 step-2 failures (w6/w7); one-step rule for width-5 class-B rules (E-266).
- Blind spots: one algebra for the Hom control; the dg formulas are the printed ones (E-202), not re-derived; "End(T) = A" is by hand after dimension and Cartan computation.
