# Scholar, conference round 108 (cycle 001, phase 2, research)

author: scholar · round: 108 · kind: conference · no new work, no scripts run

## Most promising question for phase 3
Is a K0-unimodular rigid exchange path (E-287: A_n reached from B_n at depth 3/4, n = 5, 6) a tilting mutation sequence, and does it settle T10? Two literature-checkable steps: (a) Bruning/Keller-Vossieck-style generation of tilting complexes by iterated exchange (Aihara-Iyama silting order; for a hereditary D_n piecewise-hereditary B_n the Happel-Unger / Riedtmann-Schofield exchange graph of tilting objects is connected), so "hit" = derived equivalence is a known theorem, not our search; (b) the intermediate End(T): are they LNAs? If the path leaves LNAs, the gate/guard question (H-015) reduces to a known statement: the E-187 moves are a proper subset of tilting mutations. Why this one: it is the only result that gives a certificate for A_n ~ B_n independent of the move set, it answers the "needs a non-move-set certificate" line in T10, and the sources are checkable (I own the literature part, theorist the intermediate End).

## One question outside this phase's main thread
H-025 / Phi periodic vs FCY (dormant since the plan): 1310.1557 Thm 2(a) gives FCY => Phi periodic only; no source for the converse. A single cheap test is a known piecewise hereditary Nakayama (Coxeter polynomial with non-cyclotomic factor? no: cyclotomic) with Phi periodic but not FCY. If a published example exists, H-025 is a literature citation, not a project claim. I have not searched for it since phase 1.

## Weakest claim the workshop relies on
The unsourced "derived equivalent => Phi-equivariant Z-congruence" step behind every Cartan-level NO. The congruence half is Ladkani math/0610685 Cor 3.13 / Lemma 3.14; the equivariance P Phi_A = Phi_B P is only inside the proof of Lemma 3.14, the statement is ours (R101). Nikulin Cor 1.13.3 is still UNVERIFIED (secondary only) and E-269's NO power at lam < 1.5 is untested (no certified-inequivalent pair). Runner-up: E-287's "generation (Bruning) from memory, UNVERIFIED".

## Closest to a finding this phase
E-287 (H-027, A_n reached from B_n by 3/4 tilting exchanges, End(T) = A_n verified for 7/7 and 9/9, with a control and a depth-1 reproduction of E-282). Together with Fosse Cor. 4.6 (already derived equivalence) it is consistent, but "finding" is premature. One check it still needs: confirm each exchange is a silting/tilting mutation in the Aihara-Iyama sense (one summand replaced via an exchange triangle with a minimal approximation), and that End = A_n is a full isomorphism of quivers with relations, not only thin-triangular Cartan match. The first half is a citation (Aihara-Iyama Thm 2.31/2.35) plus a triangle check on one hit; second half is toolsmith. (E-285, my own, is a finding of a different kind: numbers match; closed.)

## What I need
- theorist: whether the intermediate End(T) along E-287's depth-4 path are LNAs, and one hit written as explicit exchange triangles, so I can match each step to the Aihara-Iyama exchange-triangle statement.
- toolsmith: full quiver-with-relations isomorphism for End(T) = A_n (n = 5, 6), not Cartan match only.
- chair/human: PDFs (not on arXiv) for Happel-Seidel, BGS, Nikulin 1979, Anick 1982, Forre 2011; otherwise those citations stay secondary. Owed from me: literature files for math/0702060, 2606.15595, 2604.25789.
