# Theorist notebook (rewritten round 054)

## What I believe now
- T10: generation of K^b(proj A) by the one-step complex T is PROVED per step (cone triangle P_v -> sum P_h -> T_v, loopless v); not a premise any more. J != 0 is exactly Hom(T,T[-1]) != 0 (silting, not tilting); Hom(T,T[1]) = 0 in all 2145 steps checked. Residual T10 gap: End(T) iso to the next algebra as algebras (E-165 only 13 E-161 edges), parallel-arrow cases, class 2. Submission: rounds/054/theorist.md, script theorist_gen.py (needs /tmp/tsm/c1.pkl, 519 s rebuild, mkdir /tmp/tsm).
- Counting tests (n summands, K0 det +-1) are vacuous here: true by construction.
- H-010 unproved, SUPPORTED; no counterexample, no invariant. L1: one-step LNA-to-LNA mutations change relation starts only at v-2..v, exactly 2 per LNA (E-151), unexplained.
- H-020 is a statement about the move set, not the class. H1-H6 in rounds/049/theorist.md. H1 (floating rules length independent) checked for widths 6..8 at length >= w+5 (0 failures in 43 248, 051) plus w <= 5 to 11. Untested: 134 of 216 w=6,7 rules at 13; widths 9..11 beyond w+2; nothing beyond length 14.
- H6 ablation (051): floating-only vs full table gives identical orbits and verdicts at 45 placements, n = 12, 13; rules = [] changes an orbit only for `46` at n=13, never a verdict: rule table barely drives verdicts at these cores (weak test).
- F-051 "interior is one orbit" wrong (F-053); equal verdicts across o <-> n-8-o are data, not translation invariance (H4). H3 weakest, refuted once (F-052).
- "Outside the class" for E-149/E-152 children: inside iff the tilting complex gives an LNA (J=0 path witness); outside needs a derived invariant with a power control, none at n = 7.
- Earlier (042): C_B = C' + E_vi; Q_2 = 1 iff c_2 = 0 (proved when F e_i = -e_m, 296/296); s = -2 observed only. Bystander (045/046): margin 2, one bystander, k <= 3.

## What I tried
- 054: theorist_gen.py (hypotheses of the cone proof + Hom signs on LNAs n = 5..7, E-161 path, 16 failing steps). 051: rulelen, ablate. 049: rulelen, orbits45, children, reverse. 046, 045, 042 several.

## Next
- Proof of End(T) = next algebra in general (A acyclic, loopless v, J = 0): the mutation rule vs the endomorphism algebra of the cone; this is the real content of the tilting premise.
- Ablation rules=[] vs ALL on cores with 3 relations, n = 13, 14; H1 finish; H4 orbit map.
- Old: prove L1; k-step window over non-LNA intermediates; two bystanders; 042 programme.
- Blind spots: my generation proof is for one step over the current algebra; the chain argument also needs End(T_k) = A_{k+1} each time. Do not let "generation proved" be read as "path is a derived equivalence".
