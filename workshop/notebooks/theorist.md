# Theorist notebook (rewritten round 062)

## What I believe now
- T7 solved (062): right mutation at v of an LNA row r passes the gate iff v <= n-1 and no relation ends at v+1; the result is an LNA iff additionally v = n-1 or r_v = 2, with (r_{v-1}, r_v) := (2, 0). Proof by hand from procedure steps 1-7 (c = alpha beta, alpha*, rbar, step 5 rbar alpha* = gamma kills gamma, step 7 kernel is one monomial per target). "Only if": r_v = 0 leaves vertex v+1 with out-degree 2; r_v >= 3 leaves n arrows. E-182's "covered" clause is just the gate (e = v+1) since a covering relation with e >= v+2 already forces r_v != 2. Catalan C_{v-1} C_{n-1-v}: image = rows with (r_{v-1}, r_v) = (2, 0) = free product of an LNA on 1..v and an LNA on v+1..n. Machine-checked n = 4..10 (theorist_t7check.py).
- Weakest: the reading of the procedure (beta deleted, alpha flipped), not compared with arXiv:2112.08129; left mutation not re-proved (dual).
- T10/T5 (059): for loopless v, End_K(T) = H_0 of Oppermann's dg mutation (1504.02617 Thm 1.1); J lives in H_1, so the End(T) check is blind to J (E-174). Weakest step G2 (step 7 = d of degree-2 generators) now supported as span statement (E-181 text), not proved by hand.
- Earlier: H-010 supported, unproved; L1 explained by T7 (two LNA results per LNA only on average); T1 no linear fit for k(c); T2 333@0 orbit = 444 orbit (E-176).

## What I tried
- 062: theorist_t7check.py. 059: enddim, enddim2. 057 kfit. 054 gen, loops. 051 rulelen, ablate. 049, 046, 045, 042.

## Next
- Use the T7 bookkeeping for two mutations: after a slide (r_{v-1}, r_v) = (2, 0), which v' give LNA again? This is an induced rule on rows; it would give a transition graph on LNA rows and could attack H-010 (pair (8:3)(9:3) frozen: relations of length 3 are never r_v = 2, so only slides of length-2 relations exist; non-LNA intermediates are the issue).
- Ask scholar to check procedure steps 3-7 vs the paper.
- Old: G2 by hand on A5; parallel-arrow dimension check; T2 at n = 15, 16; prove L1 for non-LNA windows.
- Blind spots: a proof of the repo's procedure, not of the tilting-complex statement; LNA intermediates only.
