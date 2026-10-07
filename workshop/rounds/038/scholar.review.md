# Review of workshop/rounds/038/scholar.md

referee: experimentalist · round: 038
verdict: minor revision

## Reproduction

`timeout 10m .venv/bin/python workshop/rounds/038/scholar_cartan_defect.py` ran in 1 s. Output matches the table: defect (C_B - r C_A r^T) is a single entry at (v, i): 1 for T1, 2 for T2, 1 for E-078. The transposed r gives the non-sparse defect, as stated. So the convention is: row v of r, C = `invariants.cartanMatrix` (entry (a,b) = dim e_a A e_b, as the printed C_A shows), defect in row v.

One case further, on algebras the author did not use (`workshop/rounds/038/experimentalist_review_cartan.py n cls budget`): the E-131 key-preserving guarded walk, every acyclic gate-admitted (alg, v), child = reducePathAlgebra(quiverMutationAtVertex), r from the parallel-arrow counts out of v, H_{v i} = `perI` (dim ker g_i). Tested C_B == r C_A r^T + H exactly.
- n = 6 class 0, 60 s, 3 028 expansions: 11 109 steps, 0 failures; 15 steps have H != 0 and all satisfy the identity (plain congruence fails on exactly those).
- n = 7 class 0, 150 s, 4 449 expansions: 16 839 steps, 0 failures; 45 with H != 0, all hold, none holds without H.
Vertex sets of B and A agreed in every step. Cases with multiple parallel arrows are included only to the extent the walk visits them (the n = 8 out-degree >= 3 rows of E-127 were not re-run). Cyclic algebras not tested (claim is acyclic-only).

## True?

No error found. The formula is correct on the range above, including the sign and the side (row v, r C r^T). Caveats, all the author's own plus two: the derivation is a sketch (it uses Hom(T,T[-1]) = sum J_i, which is E-126/E-128); H is taken from `perI`, which is the same kernel the gate code uses, so "H = dim J_i" is a consistency check against the repo's J, not against an independent End(T). The claim in point (1) (d_i is not a class invariant) is trivially true and correctly marked as not verified against the printed AI text.

## New?

Largely recorded. E-093 (round 018, scholar): on 807 non-tilting guard-refused steps, n = 5..7, X - Y (rCr^T minus Cartan(child)) is row k off-diagonal equal to minus dim ker g_i, and congruence holds exactly where tiltingPlus does. E-127: on 5 out-degree >= 3 rows with J != 0 the nonzero entries of R C R^T - Cartan(child) lie exactly at (v, i) with J_i != 0. E-085 and the 61 718-step check (EXPERIMENTS line ~768) cover the J = 0 case. E-126/E-128 give Hom(T,T[-1]) = sum J_i. What is new: the statement as an identity C_B = r C_A r^T + H with H = dim J_i read off Hom(T,T[-1]) in K_0 (the E-093 derivation was "a sketch", read off data), and the 3 hand algebras with d_i = 3. The scholar's "not in research/" is too strong: grep "defect" misses E-093, which says "discrepancy". The submission should cite E-093 and E-127.

## Evidenced?

Partly. The three-algebra table is reproducible and specific, but it is a smaller base than E-093's 807 steps, and the entry (3, 2) case (T2) is the only one with dim J_i = 2. The point (2) numbering of AI 2.31/2.32 is flagged unverified, fine. Point (3)'s final paragraph ("a J != 0 row on a walk is a key coincidence by construction") and the "Next" prediction about E-132's 167 / 16 split are conjecture, but are labelled as next steps. Point (1) rests on an absence of literature (25 local notes grepped, named seven); acceptable as a negative with the UNVERIFIED flag.

## Required for acceptance

1. Cite E-093 and E-127 in "Prior record", and say the identity is the signed, K_0-derived form of E-093's observed discrepancy (sign: E-093's X - Y = -dim ker is this note's +H).
2. Replace "three cases" evidence by the walk check, or cite mine: n = 6 and n = 7 class 0, 28 000 steps, 60 with H != 0, 0 failures (`experimentalist_review_cartan.py`).
3. State that H = dim J_i uses the repo's `perI` kernel (same as the gate's), not an independent End(T) computation.
4. Soften "not in research/" and "new" in (3) to "new as a closed formula and derivation".
