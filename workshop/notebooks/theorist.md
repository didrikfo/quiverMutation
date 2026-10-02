# Theorist notebook (rewritten round 029)

## What I believe now
- Setup (E-110): J_i != 0 iff circuit in Gamma_i (monomial + two-term, scalar 1). W = length-2 ground path. D (nn), G, H are hand-built non-LNA-derived: their Coxeter keys are in no LNA key set (round 029, checked), so they cannot sit on a walk. That is example-specific, not a theorem.
- Round 029: skeptic's 22 n=8 c0 rows are NOT ground paths (J = 0 so no circuit): "half-W" (20/28 at (row,i) level: p1b1=p2b1, only p2 killed by b2) and "two loose pendants" (6/28). 2 involve a coefficient-2 relation (scalar != 1), unexplained.
- On walks (n=8 c0 6 776 (row,i) pairs, c1 8 990) every component of Gamma_i with >= 2 edges has exactly 2 edges; no nn, no circuit >= 3. dim e_iAe_v: {0,1,2,3} = 1020/5492/261/3 at c0; max 2 at c1. Circuit >= 3 needs dim >= 3; nn needs dim >= 2. Reduction only, not proof.
- Older: W (026) 0 mismatches on 17 802 rows; D1 depth-1 rule (022); counts are load-dependent, report by algebra count.

## What I tried
- 029: `rounds/029/theorist_circuit.py` (graph per (row,i), components with ground merged), `theorist_keys.py`. First tried the full walk at 200 s with two jobs in parallel: too slow (>6 min), use 150 s one at a time. Do not `pkill -f` with the script name inside the same command (kills the shell).
- 026: shape tables insufficient; sides-through-v key; exact match for W.

## Next
- Hom-dimension route: dim e_iAe_v = dim Hom(T_i,T_v); thin Hom for linear Nakayama modules (<= 1); find how mutation raises it and whether depth bounds it. That would prove "no circuit >= 3" (needs dim >= 3) and leave nn (dim 2).
- Is there an nn algebra whose Coxeter key is an LNA key? If none, key obstruction becomes a proof candidate.
- Prove J != 0 => W at step 7 is moot (scholar 027: false in general); instead state W as "J != 0 on LNA-derived algebras".
- Blind spots: sample is capped BFS at n=8 classes 0,1 only; scalar != 1 pairs (4) are where a counterexample could hide; "ground as one vertex" only bounds sizes.
