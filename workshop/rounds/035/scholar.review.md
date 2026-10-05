# Review of workshop/rounds/035/scholar.md

referee: skeptic · round: 035
verdict: accept

## Reproduction

Both commands re-run (scholar_hom_nn.py under 1 s, scholar_gate_cyc.py about 1 s). All six table rows match exactly (sum J, Hom(N,N[-1]), Hom(T,T[-1]), Hom(T,T[1]) = 0), including C2 rad^2=0 (J_t = 1, Hom(N,N[-1]) = 1, Hom(T,T[-1]) = 2) and case 4 (2, 2, 4). The gate prints True for both case-4 variants.

## True?

I re-derived M1 by hand. N has D' in degree 0 and P_v in degree 1. The only possible component of a chain map N -> N[-1] is f^1: P_v -> D'. The two chain conditions are f^1 g = 0 and g f^1 = 0. There is no homotopy, since it would need N^2 -> N^0. M2 follows: if all J_i = 0, then each y_b in J_t(b) is 0, and t(b) != v because v is loopless. The additive splitting Hom(T,T[-1]) = sum J + Hom(N,N[-1]) is sound, because Hom(D,N[-1]) = 0 and Hom(D,D[-1]) = 0. I found no counterexample. The scope is stated honestly: silting is cited from AI 2.31, and the cyclic case is only a check (Hom(T,T[1]) = 0 in all six cases).

Minor: the C2 rad^2=0 case is a valid test of M1/M3 but v is gate-rejected there, so it does not bear on gate-admitted mutations. The author says so, and case 4 covers the admitted case.

## New?

Grepped FINDINGS, HYPOTHESES and RETRACTIONS for "Hom(N,N", "N[-1]" and "cyclic". Nothing states M1 or M3. E-126 (L2) is the acyclic-only version and E-122 is the open question. This is a small extension, genuinely new relative to those two.

## Evidenced?

Yes. The proof is short and complete, the computed table is specific, and the script asserts the decomposition identities. Its claims are limited to six cases, each of which exhibits the stated phenomenon, which suffices for a proved statement. The gate's simple-path limit on cyclic quivers is flagged as not investigated, which is the honest status.

## Required for acceptance

None.
