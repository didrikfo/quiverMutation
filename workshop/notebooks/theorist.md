# Theorist notebook (rewritten round 057)

## What I believe now
- T10: generation of K^b(proj A) by the one-step complex T is PROVED per step (cone triangle P_v -> sum P_h -> T_v, loopless v); J != 0 is Hom(T,T[-1]) != 0 (silting, not tilting); Hom(T,T[1]) = 0 in 2145 steps. Residual gap: End(T) iso next algebra as algebras (E-165 only 13 edges), parallel arrows, class 2. Script rounds/054/theorist_gen.py (needs /tmp/tsm/c1.pkl, 519 s rebuild).
- Counting tests (n summands, K0 det +-1) are vacuous here.
- H-010 unproved, SUPPORTED. L1: one-step LNA-to-LNA mutations change relation starts only at v-2..v, exactly 2 per LNA (E-151), unexplained.
- H-020 concerns the move set, not the class. H1 floating rules length independent: checked widths 6..8 at length >= w+5; untested widths 9..11, lengths > 14. H3 weakest.
- 057 (breadth): T1 (k(c) rule): leave-one-out linear fits on 17 cores, best SSE 54, <= 6/17 exact -> no cheap rule; propose closing T1 (wording: no linear fit, 11 independent cores best SSE 48, <= 3/11 exact). T2: orbit of 333@0 and 333@(n-6) IS the 444 orbit (rowset equality n = 13..16); other 333 placements disjoint (n = 14). So E-065's "upper bound" holds only for the 33y rows, not the orbit.
- Earlier (042): C_B = C' + E_vi; Q_2 = 1 iff c_2 = 0; s = -2 observed only. Bystander (045/046): margin 2, one bystander, k <= 3.

## What I tried
- 057: theorist_kfit.py, theorist_same.py. 054: theorist_gen.py. 051: rulelen, ablate. 049: rulelen, orbits45, children, reverse. 046, 045, 042.

## Next
- T2 breakdown done n=14: 444 orbit has 2871 distinct words, only ~20 are 33y; propose closing T2 as 33y-rows statement. Run n=15,16 if wanted.
- Proof of End(T) = next algebra (A acyclic, loopless v, J = 0): the real tilting premise.
- Ablation rules=[] vs ALL on 3-relation cores n = 13, 14; H1 finish; H4 orbit map.
- Old: prove L1; k-step window; two bystanders; 042 programme.
- Blind spots: generation proof is per step; a chain also needs End(T_k) = A_{k+1}. Linear-fit failure in T1 does not show no rule exists; drift-based rule untried.
