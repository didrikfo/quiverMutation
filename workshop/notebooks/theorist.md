# Theorist notebook (rewritten round 115)

## What I believe now
- T10 / H-015 NO-capable invariant (115): HH^0..2 (Hochschild cohomology dims, derived invariant by Rickard/Happel, independent of key and Cartan) is implemented (`rounds/115/theorist_hh.py`, E-relative bar complex mod 32003). YES 13/13 on E-201's steps (consistency only: the 13 are chain-joined), NO 2/2 on planted same-Cartan pair (commutative square vs ab=0, with tail, n=4,7), NO 2/56 on refused steps. BLIND on every LNA (n<=7) and on the n=9 F-010 equal-key inequivalent pair (HH^0..6 = k in degree 0). So at n<=9 no HH-dimension test can say NO to a key-kept step in the LNA-descended region.
- Cartan-level NO for F-010 pair is q_S (E-292); there is no key- and Cartan-free NO known on any key-equal inequivalent pair. R-candidate remains YES-side only.
- Older (111): M1/M2 loop models proved (E-290, E-294); degree-0 loop cancels at w=3; M4/M5 not done; non-negative grading location unknown.

## What I tried
- 115: HH code with self-tests (A_n, Kronecker (1,3,0), D4-equivalent square), controls, all LNAs n<=7, all 75 J!=0 steps of the rebuilt walk, F-010 pair. Rebuild of the /tmp checkpoint takes ~6 min (two runs of skeptic_c1_13.py).
- Weak point: HH value (1,1,1) of ab=0 square not independently cross-checked; mod p not Q; coefficients as presented.

## Next
- Ask for a skeptic cross-check of HH (second implementation or gentle-algebra formula).
- Need a walk that leaves HH-trivial region: log HH(parent), HH(child) for every J!=0 step at n>=8 with parallel arrows; first nontrivial child HH is the chance of a real NO.
- Prove or find source: tree-quiver algebras have HH^m=0, m>=1 (would make the obstruction a theorem for LNAs).
- T5: M4/M5 g-filtration; where non-negative grading enters (from 111).
- Blind spots: dg formulas printed (E-202); Lemma A for R=k^2 only; HH over GF(p) vs Q.
